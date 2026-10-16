"""Private ctypes bridge, sharing the installed REBOUND runtime and particles."""
from __future__ import annotations

import ctypes as ct
import importlib.util
import math
import threading

import rebound

_spec = importlib.util.find_spec("rebound_bridge._bridge")
if _spec is None or _spec.origin is None:
    raise ImportError("Bridge native library not found; reinstall rebound-bridge")
lib = ct.CDLL(_spec.origin)
_lock = threading.RLock()


def _prototype(name, result, *arguments):
    function = getattr(lib, name)
    function.restype = result
    function.argtypes = list(arguments)
    return function


_prototype("reb_bridge_python_rebound_version", ct.c_char_p)
_prototype("reb_bridge_python_sizeof_simulation", ct.c_size_t)
_prototype("reb_bridge_python_sizeof_particle", ct.c_size_t)
_prototype("reb_bridge_python_offset", ct.c_size_t, ct.c_char_p)
_prototype("reb_bridge_python_bind", ct.c_int, *([ct.c_void_p] * 4))


def _check_abi():
    expected = lib.reb_bridge_python_rebound_version().decode("ascii")
    if rebound.__version__ != expected:
        raise ImportError(f"REBOUND version {rebound.__version__} must match Bridge build version {expected}")
    for typ, name in ((rebound.Simulation, "simulation"), (rebound.Particle, "particle")):
        native_size = getattr(lib, "reb_bridge_python_sizeof_" + name)()
        if ct.sizeof(typ) != native_size:
            raise ImportError(f"REBOUND {name} ABI size mismatch: Python={ct.sizeof(typ)}, Bridge={native_size}")
    aliases = {"particles": "_particles", "integrator": "_integrator", "collision": "_collision", "boundary": "_boundary"}
    fields = ("t", "G", "dt", "dt_last_done", "N", "N_var", "N_active", "particles", "exact_finish_time", "integrator", "collision", "boundary", "ri_whfast", "ri_ias15", "additional_forces", "pre_timestep_modifications", "post_timestep_modifications", "heartbeat", "extras")
    for field in fields:
        pyname = aliases.get(field, field)
        if field in ("additional_forces", "pre_timestep_modifications", "post_timestep_modifications", "heartbeat"):
            pyname = "_" + field
        if getattr(rebound.Simulation, pyname).offset != lib.reb_bridge_python_offset(field.encode("ascii")):
            raise ImportError(f"REBOUND simulation ABI mismatch at field {field}.")


_check_abi()
# 共用已加载的 REBOUND，避免混用运行库和内存分配器。
_symbols = ("reb_simulation_integrate", "reb_simulation_synchronize", "reb_simulation_free", "reb_particles_transform_jacobi_to_inertial_posvel")
if lib.reb_bridge_python_bind(*(ct.cast(getattr(rebound.clibrebound, name), ct.c_void_p) for name in _symbols)):
    raise ImportError("Could not bind the loaded REBOUND library.")

_simulation = ct.POINTER(rebound.Simulation)
_prototype("reb_bridge_create", ct.c_void_p, _simulation, ct.c_double)
_prototype("reb_bridge_free", None, ct.c_void_p)
_prototype("reb_bridge_add_subsystem", ct.c_int, ct.c_void_p, ct.c_int, _simulation, ct.c_double, ct.c_int)
_prototype("reb_bridge_python_add_restored_subsystem", ct.c_int, ct.c_void_p, ct.c_int, _simulation, ct.c_double)
_prototype("reb_bridge_advance", ct.c_int, ct.c_void_p, ct.c_double)
_prototype("reb_bridge_python_advance_steps", ct.c_int, ct.c_void_p, ct.c_double,
           ct.c_ulonglong, ct.c_ulonglong)
_prototype("reb_bridge_last_error", ct.c_int)
_prototype("reb_bridge_last_error_message", ct.c_char_p)
_prototype("reb_bridge_python_cache_valid", ct.c_int, ct.c_void_p, ct.c_int)
_prototype("reb_bridge_python_restore_cache", ct.c_int, ct.c_void_p, ct.c_int, ct.c_int)
for _name in ("method", "subsystem_coupling"):
    _prototype("reb_bridge_set_" + _name, ct.c_int, ct.c_void_p, ct.c_int)
    _prototype("reb_bridge_get_" + _name, ct.c_int, ct.c_void_p)

_METHODS = {"kdk": 0, "lf4_2": 1}
_COUPLINGS = {"monopole": 0, "resolved": 1}


class BridgeError(RuntimeError):
    """A native integration or coupling operation failed."""


def _positive(value, name):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a positive finite number.")
    return value


class NativeBridge:
    """Own only Bridge's bookkeeping; keep Python component simulations alive."""

    def __init__(self, main_sim, dt, method="lf4_2", coupling="monopole"):
        self._ptr = None
        if not isinstance(main_sim, rebound.Simulation):
            raise TypeError("main_sim must be a rebound.Simulation.")
        self.dt = _positive(dt, "Outer timestep")
        # 注意要保留引用，C 还在用这些模拟。
        self.main = main_sim
        self.subsystems = []
        self._hosts = set()
        previous_dt = main_sim.dt
        with _lock:
            self._ptr = lib.reb_bridge_create(ct.byref(main_sim), self.dt)
            main_sim.dt = previous_dt
            if not self._ptr:
                self._check(-1, "Creating Bridge")
            try:
                self.method = method
                self.coupling = coupling
            except Exception:
                self.close()
                raise

    def _require_open(self):
        if not self._ptr:
            raise RuntimeError("Bridge is closed")

    def _check(self, result, operation):
        if result != 0:
            code = lib.reb_bridge_last_error()
            detail = lib.reb_bridge_last_error_message().decode("utf-8", errors="replace")
            if code == 2:
                raise MemoryError(f"{operation}: native allocation failed.")
            if code in (1, 3, 4):
                raise ValueError(f"{operation}: {detail}")
            raise BridgeError(f"{operation}: {detail} (Bridge error {code}).")

    def _check_integration(self, result):
        # Capture Bridge's error before processing any REBOUND diagnostics.
        native_error = None
        try:
            self._check(result, "Advancing Bridge")
        except (BridgeError, ValueError, MemoryError) as exc:
            native_error = exc
        if native_error is None:
            for sim in (self.main, *self.subsystems):
                sim.process_messages()
            return

        diagnostics = []
        cause = None
        for index, sim in enumerate((self.main, *self.subsystems)):
            component = "main" if index == 0 else f"subsystem {index - 1}"
            # process_messages consumes an error before raising. Continue to
            # drain the queue so later errors cannot hide its other messages.
            while True:
                try:
                    sim.process_messages()
                except (RuntimeError, RuntimeWarning) as exc:
                    diagnostics.append(f"{component}: {exc}")
                    if cause is None:
                        cause = exc
                else:
                    break
        if diagnostics:
            native_error.args = (f"{native_error} REBOUND diagnostics: {'; '.join(diagnostics)}",)
            raise native_error from cause
        raise native_error

    @property
    def method(self):
        with _lock:
            self._require_open()
            return ("kdk", "lf4_2")[lib.reb_bridge_get_method(self._ptr)]

    @method.setter
    def method(self, value):
        self._require_open()
        if value not in _METHODS:
            raise ValueError("method must be 'kdk' or 'lf4_2'.")
        with _lock:
            self._require_open()
            self._check(lib.reb_bridge_set_method(self._ptr, _METHODS[value]), "Setting method")

    @property
    def coupling(self):
        with _lock:
            self._require_open()
            return ("monopole", "resolved")[lib.reb_bridge_get_subsystem_coupling(self._ptr)]

    @coupling.setter
    def coupling(self, value):
        self._require_open()
        if value not in _COUPLINGS:
            raise ValueError("coupling must be 'monopole' or 'resolved'.")
        with _lock:
            self._require_open()
            self._check(lib.reb_bridge_set_subsystem_coupling(self._ptr, _COUPLINGS[value]), "Setting coupling")

    def add_subsystem(self, host_index, sim, dt, *, restore=False):
        self._require_open()
        if not isinstance(sim, rebound.Simulation):
            raise TypeError("subsystem must be a rebound.Simulation.")
        if sim is self.main or any(sim is existing for existing in self.subsystems):
            raise ValueError("subsystem must be distinct from the main simulation and registered subsystems")
        if not isinstance(host_index, int) or not 0 <= host_index < self.main.N:
            raise ValueError(f"host_index {host_index!r} is outside the main particle range [0, {self.main.N})")
        if host_index in self._hosts:
            raise ValueError(f"host_index {host_index} already represents a subsystem")
        dt = _positive(dt, "Subsystem timestep")
        ratio = self.dt / dt
        if not math.isfinite(ratio) or ratio < 1 or abs(ratio - round(ratio)) > 1e-12:
            raise ValueError(f"outer timestep {self.dt:.17g} must be an integer multiple of subsystem timestep {dt:.17g}")
        if sim.G != self.main.G:
            raise ValueError(f"subsystem G={sim.G:.17g} must match main G={self.main.G:.17g}")
        if not math.isclose(sim.t, self.main.t, rel_tol=1e-15, abs_tol=1e-12):
            raise ValueError(f"subsystem time {sim.t:.17g} must match main time {self.main.t:.17g} at registration")
        with _lock:
            self._require_open()
            if restore:
                # IAS15 archives preserve the adaptive next-step proposal.
                # Registration needs the configured nominal timestep, but must
                # not replace this saved numerical state during a restart.
                saved_dt = sim.dt
                result = lib.reb_bridge_python_add_restored_subsystem(self._ptr, host_index, ct.byref(sim), dt)
                sim.dt = saved_dt
            else:
                result = lib.reb_bridge_add_subsystem(self._ptr, host_index, ct.byref(sim), dt, 0)
            self._check(result, "Registering subsystem")
            self.subsystems.append(sim)
            self._hosts.add(host_index)

    def advance(self, duration):
        self._require_open()
        duration = float(duration)
        ratio = duration / self.dt
        if not math.isfinite(ratio) or ratio < 0:
            raise ValueError("duration must be finite and non-negative.")
        if abs(ratio - round(ratio)) > 64 * math.ulp(1.0) * max(1, ratio):
            raise ValueError(f"duration={duration:g} must be an integer multiple of outer timestep {self.dt:g}.")
        with _lock:
            self._require_open()
            self._check_integration(lib.reb_bridge_advance(self._ptr, duration))

    def advance_steps(self, epoch, completed, count):
        """Advance a counted batch on the original absolute outer-step grid."""
        self._require_open()
        epoch = float(epoch)
        if not math.isfinite(epoch):
            raise ValueError("Initial epoch must be finite.")
        if (type(completed) is not int or type(count) is not int
                or completed < 0 or count < 0
                or completed + count > 2**53 - 1):
            raise ValueError("completed and count must be non-negative integers with a sum <= 2**53 - 1")
        if not math.isfinite(epoch + (completed + count) * self.dt):
            raise ValueError("Resulting integration time must be finite.")
        with _lock:
            self._require_open()
            self._check_integration(lib.reb_bridge_python_advance_steps(
                self._ptr, epoch, completed, count))

    def close(self):
        with _lock:
            if self._ptr:
                lib.reb_bridge_free(self._ptr)
                self._ptr = None

    def cache_flags(self):
        """Capture valid Jacobi snapshots, with the main simulation first."""
        self._require_open()
        with _lock:
            self._require_open()
            return [bool(lib.reb_bridge_python_cache_valid(self._ptr, i)) for i in range(-1, len(self.subsystems))]

    def restore_caches(self, flags):
        """Restore saved snapshots without changing Cartesian/Jacobi values."""
        self._require_open()
        if len(flags) != 1 + len(self.subsystems) or any(type(flag) is not bool for flag in flags):
            raise ValueError("Cache flags must contain one boolean per component, main first.")
        with _lock:
            self._require_open()
            for i, flag in enumerate(flags):
                if lib.reb_bridge_python_restore_cache(self._ptr, i - 1, flag):
                    raise BridgeError(f"Cannot restore WHFast cache for component {i}")

    def __del__(self):
        if getattr(self, "_ptr", None):
            self.close()
