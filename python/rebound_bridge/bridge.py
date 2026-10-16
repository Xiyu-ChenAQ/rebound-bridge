"""A REBOUND-style interface to fixed multirate hierarchies."""

import math
from collections.abc import Mapping
from types import MappingProxyType

import rebound

from ._native import NativeBridge
from ._integration import validate_integration_request
from .checkpoint import _BRIDGE_SCHEMA_VERSION, load_checkpoint, save_checkpoint
from .hierarchy import (
    Subsystem, _ADAPTIVE_INTEGRATORS, _COORDINATES, check_simulation,
    checked_integrator, copy_settings, particle_index, positive, synchronized_copy,
)


class Bridge:
    """Couple a main REBOUND simulation to fixed, independently integrated systems.
    """

    def __init__(self, simulation=None, *, dt, method="lf4_2", coupling="monopole"):
        self._dt = positive(dt, "dt")
        self._native = None
        self._closed = False
        self._failed = False
        self._subsystems = {}
        self._host_identities = {}
        self._tree_metadata = None
        if simulation is None:
            self._main = rebound.Simulation()
            self._main.integrator = "whfast"
            self._main.dt = self.dt
        elif isinstance(simulation, rebound.Simulation):
            check_simulation(simulation)
            self._main = synchronized_copy(simulation)
        else:
            raise TypeError("simulation must be a rebound.Simulation")
        self._epoch = self._main.t
        self._completed_steps = 0
        self.method = method
        self.coupling = coupling

    @property
    def main(self):
        return self._main

    @property
    def subsystems(self):
        return MappingProxyType(self._subsystems)

    @property
    def t(self):
        return self._main.t

    @property
    def dt(self):
        return self._dt

    @property
    def method(self):
        return self._method

    @method.setter
    def method(self, value):
        if value not in ("kdk", "lf4_2"):
            raise ValueError("method must be 'kdk' or 'lf4_2'")
        if self._native is not None:
            self._native.method = value
        self._method = value

    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, value):
        if value not in ("monopole", "resolved"):
            raise ValueError("coupling must be 'monopole' or 'resolved'")
        if self._native is not None:
            self._native.coupling = value
        self._coupling = value

    def _require_open(self):
        if self._closed:
            raise RuntimeError("Bridge is closed")
        if self._failed:
            raise RuntimeError("Integration was interrupted or failed; restore the last completed checkpoint")

    def _require_editable(self):
        self._require_open()
        if self._native is not None:
            raise RuntimeError("Hierarchy membership is fixed after build; create a new Bridge to change it")

    def add(self, particle=None, **kwargs):
        """Add a main-system particle with REBOUND's usual arguments."""
        self._require_editable()
        if kwargs.get("primary") is not None and not isinstance(kwargs["primary"], rebound.Particle):
            kwargs["primary"] = self.main.particles[kwargs["primary"]]
        self.main.add(particle, **kwargs)

    def add_subsystem(self, primary=None, *, dt, name=None, integrator="whfast", **host):
        """Resolve an existing host, or add one with REBOUND particle keywords.
        """
        self._require_editable()
        integrator = checked_integrator(integrator)
        dt = positive(dt, "subsystem dt")
        ratio = self.dt / dt
        if not math.isfinite(ratio) or ratio < 1 or abs(ratio - round(ratio)) > 1e-12:
            raise ValueError(f"outer dt={self.dt:.17g} must be an integer multiple of subsystem dt={dt:.17g}")
        if primary is None:
            if not host:
                raise ValueError("primary or host particle arguments are required")
            new_host = rebound.Particle(simulation=self.main, **host)
            index = self.main.N
            default_name = host.get("hash", f"subsystem_{index}")
        else:
            if host:
                raise ValueError("primary and host particle arguments cannot be supplied together")
            index = particle_index(self.main, primary)
            new_host = None
            default_name = primary if isinstance(primary, str) else f"subsystem_{index}"
        name = str(default_name if name is None else name)
        if not name:
            raise ValueError("Subsystem name must not be empty")
        host_indices = self._host_indices()
        if name in self._subsystems or any(self._resolve_host_index(s, host_indices) == index for s in self._subsystems.values()):
            raise ValueError(f"Host or subsystem {name!r} is already registered")
        planet = (new_host if new_host is not None else self.main.particles[index]).copy()
        if not math.isfinite(planet.m) or planet.m <= 0:
            raise ValueError("host mass must be finite and positive")
        # 宿主先放局部原点，build 时再转质心系
        for attr in _COORDINATES:
            setattr(planet, attr, 0)
        sim = rebound.Simulation()
        copy_settings(self.main, sim)
        sim.integrator = integrator
        sim.dt = dt
        sim.add(planet)
        sub = Subsystem(self, name, index, sim, dt)
        if new_host is not None:
            self.main.add(new_host)
        original = self.main.particles[index]
        key = original.hash.value
        occurrences = len(host_indices.get(key, ())) + (new_host is not None)
        unique_key = key if key and occurrences == 1 else None
        # Store values, never a particle view into REBOUND's reallocating array.
        self._host_identities[name] = (unique_key, self._host_state(original), self.main.N)
        self._subsystems[name] = sub
        return sub

    @staticmethod
    def _host_state(particle):
        return (particle.hash.value, particle.m, *(getattr(particle, attr) for attr in _COORDINATES))

    def _host_indices(self):
        indices = {}
        for index, particle in enumerate(self.main.particles):
            indices.setdefault(particle.hash.value, []).append(index)
        return indices

    def _resolve_host_index(self, sub, host_indices):
        key, original, count = self._host_identities[sub.name]
        if key is not None:
            indices = host_indices.get(key, ())
            if len(indices) != 1:
                raise ValueError(f"Host {sub.name!r} was removed or no longer has a unique hash")
            return indices[0]
        # 没有唯一 hash 时，只接受原槽位和状态都未改变的宿主。
        if (self.main.N < count or sub.host_index >= self.main.N
                or self._host_state(self.main.particles[sub.host_index]) != original):
            raise ValueError(
                f"Cannot identify host {sub.name!r} after particle edits; "
                "assign a unique nonzero hash before registration"
            )
        return sub.host_index

    @classmethod
    def from_simulation(cls, simulation, *, subsystems, dt, dt_sub, method="lf4_2", coupling="monopole", integrator="whfast"):
        if not isinstance(subsystems, Mapping):
            raise TypeError("subsystems must map names to lists of particle indices or hashes")
        b = cls(simulation, dt=dt, method=method, coupling=coupling)
        groups, used = [], set()
        for name, members in subsystems.items():
            indices = [particle_index(b.main, member) for member in members]
            if not indices or len(set(indices)) != len(indices) or used.intersection(indices):
                raise ValueError(f"Subsystem {name!r} members must be nonempty, unique and disjoint from other subsystems")
            used.update(indices)
            groups.append((name, indices, [b.main.particles[i].copy() for i in indices]))
        removed = {i for _, indices, _ in groups for i in indices[1:]}
        for index in sorted(removed, reverse=True):
            b.main.remove(index=index, keep_sorted=True)
        for name, indices, particles in groups:
            host_index = indices[0] - sum(i < indices[0] for i in removed)
            sub_dt = dt_sub[name] if isinstance(dt_sub, Mapping) else dt_sub
            sub = b.add_subsystem(primary=host_index, name=name, dt=sub_dt, integrator=integrator)
            for member in particles[1:]:
                for attr in _COORDINATES:
                    setattr(member, attr, getattr(member, attr) - getattr(particles[0], attr))
                sub.add(member)
        return b

    def build(self):
        """Prepare COM proxies and register the fixed hierarchy once."""

        self._require_open()
        if self._native is not None:
            self._validate_topology()
            return self
        check_simulation(self.main)
        if not self.main.N:
            raise ValueError("main simulation has no particles")
        main = synchronized_copy(self.main)
        host_indices = self._host_indices()
        components = []
        for sub in self._subsystems.values():
            host_index = self._resolve_host_index(sub, host_indices)
            check_simulation(sub.simulation)
            if sub.simulation.dt != sub.dt:
                raise ValueError(
                    f"Subsystem {sub.name!r}: simulation.dt={sub.simulation.dt:.17g} differs from dt={sub.dt:.17g}; "
                    "set subsystem.dt before build"
                )
            if sub.simulation.G != main.G or sub.simulation.t != main.t or sub.simulation.units != main.units:
                raise ValueError(f"Subsystem {sub.name!r} must use the same units, G and epoch as the main simulation")
            sim = synchronized_copy(sub.simulation)
            center = sim.com()
            proxy = main.particles[host_index]
            if proxy.m != sim.particles[0].m:
                raise ValueError(f"Host mass for {sub.name!r} changed after adding its subsystem")
            # center 相对宿主，代理要用全局质心
            proxy.m = center.m
            for attr in _COORDINATES:
                setattr(proxy, attr, getattr(proxy, attr) + getattr(center, attr))
            sim.move_to_com()
            components.append((sub, sim, host_index))
        native = NativeBridge(main, self.dt, method=self.method, coupling=self.coupling)
        try:
            for sub, sim, host_index in components:
                native.add_subsystem(host_index, sim, sub.dt)
        except BaseException:
            native.close()
            raise
        self._main = main
        for sub, sim, host_index in components:
            sub._simulation = sim
            sub._host_index = host_index
        self._native = native
        self._epoch = main.t
        self._topology = self._topology_signature()
        self._configuration = self._configuration_signature()
        self._accepted_times = self._time_signature()
        return self

    def _topology_signature(self):
        return tuple(tuple((p.hash.value, p.m) for p in sim.particles) for sim in self._simulations().values())

    def _simulations(self):
        return {"main": self.main, **{f"sub_{i}": s.simulation for i, s in enumerate(self._subsystems.values())}}

    def _configuration_signature(self):
        # Adaptive dt is the solver's next-step proposal, not a fixed setting.
        return tuple((sim.integrator, sim.G, tuple(sim.units.values()),
                      None if sim.integrator in _ADAPTIVE_INTEGRATORS else sim.dt,
                      sim.ri_whfast.coordinates if sim.integrator == "whfast" else None)
                     for sim in self._simulations().values())

    def _time_signature(self):
        return tuple(sim.t for sim in self._simulations().values())

    def _validate_topology(self):
        if self._topology_signature() != self._topology:
            raise RuntimeError("Particle masses or topology changed after registration")
        if self._time_signature() != self._accepted_times:
            raise RuntimeError("Component time changed outside Bridge; use Bridge.integrate to advance the hierarchy")
        for sim in self._simulations().values():
            check_simulation(sim)
            if sim.G != self.main.G or sim.units != self.main.units or not math.isclose(sim.t, self.t, rel_tol=0, abs_tol=1e-12):
                raise ValueError("All component simulations must retain common units, G and time")
        if self._configuration_signature() != self._configuration:
            raise RuntimeError("Component configuration changed after build: integrator, coordinates, units, G or fixed dt")
        for sub in self._subsystems.values():
            if sub.integrator not in _ADAPTIVE_INTEGRATORS and sub.simulation.dt != sub.dt:
                raise ValueError(f"Subsystem {sub.name!r} dt disagrees with its configured fixed timestep")
        expected = self._epoch + self._completed_steps * self.dt
        if not math.isclose(self.t, expected, rel_tol=0, abs_tol=1e-11 * max(1, abs(expected))):
            raise RuntimeError("Component time changed outside Bridge; use Bridge.integrate to advance the hierarchy")

    def _accept_native_mass_roundoff(self):
        current = self._topology_signature()
        if current == self._topology:
            return
        for sim, before, after in zip(self._simulations().values(), self._topology, current):
            if before == after:
                continue
            total_mass = sum(mass for _, mass in before)
            # The total and subsequent subtractions take at most 2*N rounded
            # operations. Also cap the allowance at the central-mass scale:
            # cancellation with a much larger companion is not harmless drift.
            bound = (8 * len(before) * min(math.ulp(total_mass), math.ulp(before[0][1]))) if before else 0
            if not (sim.integrator == "whfast" and sim.ri_whfast.coordinates == "barycentric"
                    and len(before) == len(after) and before and before[1:] == after[1:]
                    and before[0][0] == after[0][0] and math.isfinite(total_mass)
                    and math.isfinite(after[0][1]) and after[0][1] > 0
                    and abs(after[0][1] - before[0][1]) <= bound):
                raise RuntimeError(
                    "Particle masses or topology changed beyond roundoff during integration; "
                    "check particle edits and WHFast barycentric mass reconstruction"
                )
        self._topology = current

    def _target_steps(self, t):
        return validate_integration_request(
            t, epoch=self._epoch if self._native is not None else self.main.t,
            dt=self.dt, completed_steps=self._completed_steps,
        )[0]

    def integrate(self, t, *, checkpoint=None, checkpoint_interval=None):
        """Advance to an absolute time on the outer-step grid.

        Optionally replace a checkpoint after every ``checkpoint_interval``
        of simulated time. A crash or interruption leaves the last successfully
        written file available; load it with ``Bridge.from_file``.
        """
        self._require_open()
        if self._native is not None and float(t) == self.t:
            # A supported legacy checkpoint may retain historical clock
            # roundoff. Requesting that exact current time is still a no-op;
            # the trusted-time check below rejects external clock edits.
            t = self._epoch + self._completed_steps * self.dt
        # 用整数步数计时，避免反复累加 t
        target, cadence = validate_integration_request(
            t, epoch=self._epoch if self._native is not None else self.main.t,
            dt=self.dt, completed_steps=self._completed_steps,
            checkpoint=checkpoint, checkpoint_interval=checkpoint_interval,
        )
        self.build()
        self._validate_topology()
        if checkpoint is not None:
            self.save_to_file(checkpoint)
        while self._completed_steps < target:
            count = min(cadence, target - self._completed_steps)
            self._validate_topology()
            try:
                self._native.advance_steps(self._epoch, self._completed_steps, count)
                self._accept_native_mass_roundoff()
                self._accepted_times = self._time_signature()
                self._completed_steps += count
                # C 正常返回后仍需检查 NaN 和溢出。
                self._validate_topology()
            except BaseException:
                # 可能已推进一部分，失败后从检查点恢复
                self._failed = True
                raise
            if checkpoint is not None:
                self.save_to_file(checkpoint)

    def to_rebound(self):
        """Return a detached simulation of physical particles in global coordinates."""
        self.build()
        result = rebound.Simulation()
        copy_settings(self.main, result)
        by_host = {s.host_index: s for s in self._subsystems.values()}
        # 加回质心位置和速度，代理不重复导出
        for i, particle in enumerate(self.main.particles):
            if i not in by_host:
                result.add(particle)
            else:
                for local in by_host[i].particles:
                    physical = local.copy()
                    for attr in _COORDINATES:
                        setattr(physical, attr, getattr(local, attr) + getattr(particle, attr))
                    result.add(physical)
        return result

    def save_to_file(self, filename):
        """Atomically save the whole hierarchy at a completed Bridge boundary."""
        self.build()
        self._validate_topology()
        metadata = {
            "bridge_schema": _BRIDGE_SCHEMA_VERSION, "dt": self.dt, "method": self.method, "coupling": self.coupling,
            "epoch": self._epoch, "completed_steps": self._completed_steps,
            "subsystems": [{"name": s.name, "host_index": s.host_index, "dt": s.dt} for s in self._subsystems.values()],
            "cache_flags": self._native.cache_flags(),
        }
        if self._tree_metadata is not None:
            metadata["physical_tree"] = self._tree_metadata
        save_checkpoint(filename, self._simulations(), metadata)

    @classmethod
    def from_file(cls, filename):
        """Restore a checkpoint made by ``save_to_file`` or autosaving integration."""
        from .checkpoint import CheckpointError
        try:
            sims, metadata = load_checkpoint(filename, require_bridge_schema=True)
            records = metadata["subsystems"]
            if not isinstance(records, list) or set(sims) != {"main", *(f"sub_{i}" for i in range(len(records)))}:
                raise ValueError("Checkpoint components do not match the hierarchy")
            # 按存档布局恢复，不重复做质心变换
            b = cls.__new__(cls)
            b._native = None
            b._closed = False
            b._failed = False
            b._dt = positive(metadata["dt"], "dt")
            b._main = sims["main"]
            b._subsystems = {}
            b._host_identities = {}
            b._tree_metadata = metadata.get("physical_tree")
            b._epoch = float(metadata["epoch"])
            b._completed_steps = metadata["completed_steps"]
            if not math.isfinite(b._epoch) or type(b._completed_steps) is not int or b._completed_steps < 0:
                raise ValueError("Invalid checkpoint time or step count")
            expected = b._epoch + b._completed_steps * b.dt
            if not math.isfinite(expected):
                raise ValueError("Invalid checkpoint time or step count")
            if abs(b.t - expected) > b.dt / 32:
                raise ValueError(
                    "Checkpoint time drift exceeds dt/32 from the recorded outer-step grid; "
                    "inspect the file with the version that created it"
                )
            b.method = metadata["method"]
            b.coupling = metadata["coupling"]
            indices = set()
            for i, record in enumerate(records):
                name, index = record["name"], record["host_index"]
                if not isinstance(name, str) or not name or name in b._subsystems or type(index) is not int or index in indices or not 0 <= index < b.main.N:
                    raise ValueError("Invalid or duplicate subsystem identity")
                indices.add(index)
                b._subsystems[name] = Subsystem(b, name, index, sims[f"sub_{i}"], positive(record["dt"], "subsystem dt"))
            b._topology = b._topology_signature()
            b._configuration = b._configuration_signature()
            b._accepted_times = b._time_signature()
            b._validate_topology()
            b._native = NativeBridge(b.main, b.dt, method=b.method, coupling=b.coupling)
            for sub in b._subsystems.values():
                b._native.add_subsystem(sub.host_index, sub.simulation, sub.dt, restore=True)
            b._native.restore_caches(metadata["cache_flags"])
            return b
        except Exception as exc:
            if "b" in locals() and b._native is not None:
                b._native.close()
            raise CheckpointError(f"Invalid Bridge checkpoint: {exc}") from exc

    def close(self):
        if self._native is not None:
            self._native.close()
        self._closed = True

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()

    def __repr__(self):
        return f"Bridge(t={self.t:g}, dt={self.dt:g}, method={self.method!r}, subsystems={list(self._subsystems)!r})"
