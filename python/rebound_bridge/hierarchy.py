"""Prepare fixed host/member systems using ordinary REBOUND particles."""

import math

import rebound

from .checkpoint import _check_serializable


_COORDINATES = ("x", "y", "z", "vx", "vy", "vz")
_ORBIT_ARGUMENTS = {
    "a", "P", "e", "inc", "Omega", "omega", "pomega", "f", "M", "E",
    "l", "theta", "T", "h", "k", "ix", "iy", "pal_h", "pal_k", "pal_ix", "pal_iy",
}
_SUPPORTED_INTEGRATORS = frozenset({"whfast", "ias15", "bs", "leapfrog", "saba", "eos"})
_ADAPTIVE_INTEGRATORS = frozenset({"ias15", "bs"})


def checked_integrator(value):
    """Normalize canonical names without silently discarding alias settings."""
    if isinstance(value, str) and value.lower() in _SUPPORTED_INTEGRATORS:
        return value.lower()
    probe = rebound.Simulation()
    probe.integrator = value
    name = probe.integrator
    if name not in _SUPPORTED_INTEGRATORS:
        raise ValueError(
            f"Integrator {name!r} is not supported; "
            "use 'whfast', 'ias15', 'bs', 'leapfrog', 'saba', or 'eos'"
        )
    if isinstance(value, str) and value.lower() != name:
        raise ValueError(
            f"Integrator alias {value!r} is not supported; use canonical name {name!r} "
            "and configure options on Bridge.main or Subsystem.simulation"
        )
    return name


def synchronized_copy(source):
    """Copy a component and expose its physical state without touching the input."""
    result = source.copy()
    if result.integrator == "whfast":
        keep, dt = result.ri_whfast.keep_unsynchronized, result.dt
        try:
            # Commit the pending half drift/corrector in the private copy.
            # It belongs to the last executed step, not a newly configured dt.
            result.ri_whfast.keep_unsynchronized = 0
            if not result.ri_whfast.is_synchronized and result.dt_last_done != 0:
                result.dt = result.dt_last_done
            result.synchronize()
        finally:
            result.dt = dt
            result.ri_whfast.keep_unsynchronized = keep
    else:
        result.synchronize()
    return result


def positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}")
    return value


def check_simulation(sim):
    """Reject states whose meaning would be lost when copying or partitioning."""
    # 只接受已适配的配置
    checked_integrator(sim.integrator)
    if sim.integrator == "whfast" and sim.ri_whfast.kernel not in ("default", "composition"):
        raise ValueError("WHFast ri_whfast.kernel must be 'default' or 'composition'")
    if sim.integrator == "saba" and (sim.ri_saba.safe_mode != 1 or sim.ri_saba.keep_unsynchronized != 0):
        raise ValueError("SABA requires ri_saba.safe_mode=1 and ri_saba.keep_unsynchronized=0")
    if sim.integrator == "saba" and sim.ri_saba._type not in range(10):
        raise ValueError("Only uncorrected SABA compositions are supported")
    if sim.integrator == "eos" and sim.ri_eos.safe_mode != 1:
        raise ValueError("EOS requires ri_eos.safe_mode=1")
    positive(sim.G, "G")
    positive(sim.dt, "component dt")
    if sim.exact_finish_time != 1:
        raise ValueError("Component exact_finish_time must be 1")
    if not math.isfinite(sim.t):
        raise ValueError("Simulation time must be finite")
    # Fixed components use half-step times in the native segment scheduler
    # (and WHFast/leapfrog advance their own clocks by half steps). Python's
    # autonomous adaptive components use a local clock within each stage.
    if sim.integrator not in _ADAPTIVE_INTEGRATORS and sim.t + sim.dt / 2 == sim.t:
        raise ValueError(
            f"Component dt/2 is below time resolution at t={sim.t!r}; "
            "use a relative initial epoch"
        )
    if sim.N_var or sim.N_active not in (-1, sim.N):
        raise ValueError("N_var must be 0 and N_active must be -1 or N")
    if sim.collision != "none" or sim.boundary != "none" or sim.softening != 0:
        raise ValueError("collision and boundary must be 'none'; softening must be 0")
    if sim.gravity not in ("basic", "compensated"):
        raise ValueError("gravity must be 'basic' or 'compensated'")
    if sim._calculate_megno:
        raise ValueError("MEGNO is not supported")
    if sim._display_data or sim._server_data:
        raise ValueError("Live REBOUND display or server state is not supported")
    if (sim._simulationarchive_filename or sim.simulationarchive_auto_interval
            or sim.simulationarchive_auto_walltime or sim.simulationarchive_auto_step):
        raise ValueError("Native Simulationarchive output is not supported; use Bridge.save_to_file()")
    _check_serializable("simulation", sim)
    for i, p in enumerate(sim.particles):
        if not math.isfinite(p.m) or p.m < 0:
            raise ValueError(f"Particle {i} mass must be finite and non-negative")
        if not all(math.isfinite(getattr(p, attr)) for attr in _COORDINATES):
            raise ValueError(f"Particle {i} position and velocity must be finite")
    # WHFast's coordinate transforms (also used by SABA) divide by the
    # first particle's mass. EOS uses Cartesian drifts and does not need this.
    if sim.integrator in ("whfast", "saba") and sim.N and sim.particles[0].m == 0:
        raise ValueError(
            f"{sim.integrator.upper()} requires the first particle to have positive mass; "
            "add the central massive body first or use IAS15"
        )


def particle_index(sim, selector):
    if isinstance(selector, rebound.Particle):
        # A hash gives stable identity even if the particle array was reallocated.
        if selector.hash.value:
            selector = selector.hash
        else:
            index = selector.index
            if index < 0:
                raise ValueError("Use a particle index or hash to identify the host")
            selector = index
    try:
        return sim.particles[selector].index
    except (KeyError, IndexError, rebound.ParticleNotFound, TypeError) as exc:
        raise ValueError(f"No particle matches {selector!r}") from exc


def copy_settings(source, target):
    # 先设单位，再恢复 G；积分器和步长另设
    if all(source.units.values()):
        target.units = tuple(source.units.values())
    target.G = source.G
    target.t = source.t


class Subsystem:
    """A host and its members, stored in local coordinates.

    Before ``Bridge.build()``, the host is at the local origin. Afterwards
    ``simulation`` uses its centre-of-mass frame. Use ``Bridge.to_rebound()``
    for global physical coordinates.
    """

    def __init__(self, owner, name, host_index, simulation, dt):
        # 保留 owner，编辑时检查 Bridge 状态
        self._owner = owner
        self._name = name
        self._host_index = host_index
        self._simulation = simulation
        self._dt = dt

    @property
    def name(self):
        return self._name

    @property
    def host_index(self):
        return self._host_index

    @property
    def simulation(self):
        return self._simulation

    @property
    def dt(self):
        return self._dt

    @dt.setter
    def dt(self, value):
        self._owner._require_editable()
        value = positive(value, "subsystem dt")
        # 修改步长也要检查整数比，留一点浮点容差
        ratio = self._owner.dt / value
        if not math.isfinite(ratio) or ratio < 1 or abs(ratio - round(ratio)) > 1e-12:
            raise ValueError(f"Outer dt / subsystem dt must be a positive integer; got {ratio!r}")
        self._dt = value
        self.simulation.dt = value

    @property
    def particles(self):
        return self.simulation.particles

    @property
    def integrator(self):
        return self.simulation.integrator

    @integrator.setter
    def integrator(self, value):
        self._owner._require_editable()
        self.simulation.integrator = checked_integrator(value)

    def add(self, particle=None, **kwargs):
        """Add a local particle, defaulting orbital elements to the host primary."""
        self._owner._require_editable()
        if isinstance(particle, str) or (isinstance(particle, list) and any(isinstance(p, str) for p in particle)):
            raise ValueError("Horizons states are global; add them to rebound.Simulation and use Bridge.from_simulation()")
        # 默认绕宿主，不绕已有成员的质心
        primary = kwargs.get("primary")
        if primary is not None and not isinstance(primary, rebound.Particle):
            kwargs["primary"] = self.particles[primary]
        elif primary is None and any(kwargs.get(key) is not None for key in _ORBIT_ARGUMENTS):
            kwargs["primary"] = self.particles[0]
        self.simulation.add(particle, **kwargs)

    def __repr__(self):
        return f"Subsystem(name={self.name!r}, N={self.simulation.N}, dt={self.dt:g})"
