"""Atomic Bridge checkpoints containing JSON and native REBOUND archives."""

from collections.abc import Mapping
import ctypes
from functools import lru_cache
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import sys
import tempfile
import warnings
import zipfile

import rebound


# 积分器内部状态一起存，只有位置速度不够续算。
_FORMAT = "rebound_bridge_checkpoint"
_FORMAT_VERSION = 1
# These describe on-disk semantics, independently of the package release.
# Bump the relevant schema and add an explicit reader/migration whenever those
# semantics change, including the meaning of saved integration/cache state.
_BRIDGE_SCHEMA_VERSION = 1
_TREE_SCHEMA_VERSION = 1
_MAX_MANIFEST_BYTES = 8 * 1024 * 1024


class CheckpointError(ValueError):
    """A checkpoint is invalid or incompatible with the current hierarchy."""


def _abi():
    return {
        "system": platform.system(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "pointer_size": ctypes.sizeof(ctypes.c_void_p),
        "long_size": ctypes.sizeof(ctypes.c_long),
        "simulation_size": ctypes.sizeof(rebound.Simulation),
        "particle_size": ctypes.sizeof(rebound.Particle),
    }


@lru_cache(maxsize=1)
def _bs_nbody_derivative_address():
    """Identify REBOUND's private native callback without trusting user state."""
    # nbody_derivatives is a static C function, so it cannot be resolved through
    # clibrebound. One stationary particle creates its default ODE cheaply. Only
    # the native function address is retained; the temporary simulation is freed.
    probe = rebound.Simulation()
    probe.integrator = "bs"
    probe.dt = 1e-6
    probe.add(m=1)
    probe.step()
    return ctypes.cast(probe.ri_bs._nbody_ode.contents._derivatives, ctypes.c_void_p).value


def _has_external_odes(sim):
    if not sim._N_odes:
        return False
    # BS uses the public ODE registry for its own N-body work array. Native
    # archives rebuild this array from the saved particles on the next step.
    # User ODEs and changes to the default callback/model are not reconstructible.
    nbody = sim.ri_bs._nbody_ode
    if sim.integrator != "bs" or sim._N_odes != 1 or not sim._odes or not nbody:
        return True
    if ctypes.cast(sim._odes[0], ctypes.c_void_p).value != ctypes.cast(nbody, ctypes.c_void_p).value:
        return True
    ode = nbody.contents
    return (
        ode.length != 6 * sim.N
        or bool(ode.needs_nbody)
        or bool(ode.ref)
        or ctypes.cast(ode.r, ctypes.c_void_p).value != ctypes.addressof(sim)
        or bool(ode._getscale)
        or bool(ode._pre_timestep)
        or bool(ode._post_timestep)
        or ctypes.cast(ode._derivatives, ctypes.c_void_p).value != _bs_nbody_derivative_address()
    )


def _check_serializable(name, sim):
    if not isinstance(sim, rebound.Simulation):
        raise TypeError(f"Checkpoint component {name!r} is not a rebound.Simulation.")
    # callback 只有地址，不能随存档恢复。
    callbacks = [
        field.lstrip("_")
        for field, field_type in sim._fields_
        if issubclass(field_type, ctypes._CFuncPtr) and bool(getattr(sim, field))
    ]

    default_mercurius = ctypes.cast(
        rebound.clibrebound.reb_integrator_mercurius_L_mercury, ctypes.c_void_p
    ).value
    for structure_name in ("ri_mercurius", "ri_trace"):
        structure = getattr(sim, structure_name)
        for field, field_type in structure._fields_:
            if not issubclass(field_type, ctypes._CFuncPtr):
                continue
            callback = getattr(structure, field)
            if not callback:
                continue
            if structure_name == "ri_mercurius" and field == "_L":
                if ctypes.cast(callback, ctypes.c_void_p).value == default_mercurius:
                    continue
            callbacks.append(f"{structure_name}.{field.lstrip('_')}")
    if callbacks:
        raise ValueError(
            f"Cannot copy or checkpoint component {name!r}: "
            f"REBOUND cannot restore callbacks {', '.join(callbacks)}"
        )
    if sim.extras or _has_external_odes(sim) or any(p.ap for p in sim.particles):
        raise ValueError(
            f"Cannot copy or checkpoint component {name!r}: "
            "external extras, ODEs or particle.ap state cannot be restored"
        )


def _simulation_bytes(sim):
    # Follow REBOUND 4.6 Simulation.__reduce__'s native stream allocation/free
    # protocol, without invoking pickle or storing executable Python objects.
    buffer = ctypes.c_char_p()
    size = ctypes.c_size_t()
    rebound.clibrebound.reb_simulation_save_to_stream(
        ctypes.byref(sim), ctypes.byref(buffer), ctypes.byref(size)
    )
    try:
        if not buffer or not size.value:
            raise RuntimeError("REBOUND could not serialize the simulation.")
        return ctypes.string_at(buffer, size.value)
    finally:
        if buffer:
            rebound.clibrebound.reb_simulation_output_free_stream(buffer)


def _json_bytes(value):
    try:
        return json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError("Checkpoint metadata must be JSON-serializable with no NaN or infinity") from exc


def save_checkpoint(path, simulations, metadata):

    if not isinstance(simulations, Mapping) or not simulations:
        raise ValueError("simulations must be a nonempty mapping")
    if not isinstance(metadata, dict):
        raise TypeError("Checkpoint metadata must be a dictionary.")
    _json_bytes(metadata)
    for name, sim in simulations.items():
        if not isinstance(name, str) or not name:
            raise ValueError("Checkpoint component names must be nonempty strings.")
        _check_serializable(name, sim)

    entries = []
    payloads = {}
    for index, (name, sim) in enumerate(simulations.items()):
        member = f"simulations/{index}.bin"
        payload = _simulation_bytes(sim)
        payloads[member] = payload
        entries.append({
            "name": name,
            "member": member,
            "size": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        })
    manifest = _json_bytes({
        "format": _FORMAT,
        "format_version": _FORMAT_VERSION,
        "bridge_version": version("rebound-bridge"),
        "rebound_version": rebound.__version__,
        "abi": _abi(),
        "simulations": entries,
        "metadata": metadata,
    })
    if len(manifest) > _MAX_MANIFEST_BYTES:
        raise ValueError("Checkpoint metadata is too large.")

    target = Path(path)
    temporary = None
    # 同目录写临时文件，完整写入后再替换旧存档。
    try:
        with tempfile.NamedTemporaryFile(
            mode="w+b", prefix=f".{target.name}.", suffix=".tmp",
            dir=target.parent, delete=False,
        ) as output:
            temporary = Path(output.name)
            with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
                archive.writestr("manifest.json", manifest)
                for member, payload in payloads.items():
                    archive.writestr(member, payload)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, target)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _reject_nonfinite(value):
    raise ValueError(f"Non-finite JSON value in checkpoint: {value}")


def _check_schema(value, expected, label):
    if type(value) is not int or value != expected:
        raise ValueError(
            f"Unsupported {label} {value!r}; supported version is {expected}"
        )


def load_checkpoint(path, *, require_bridge_schema=False):

    payloads = []
    try:
        with zipfile.ZipFile(path, "r") as archive:
            names = archive.namelist()
            if len(names) != len(set(names)):
                raise ValueError("Checkpoint contains duplicate archive members.")
            if archive.getinfo("manifest.json").file_size > _MAX_MANIFEST_BYTES:
                raise ValueError("Checkpoint manifest is too large.")
            manifest = json.loads(archive.read("manifest.json"), parse_constant=_reject_nonfinite)
            if not isinstance(manifest, dict) or manifest.get("format") != _FORMAT:
                raise ValueError("This is not a REBOUND_BRIDGE checkpoint.")
            _check_schema(manifest.get("format_version"), _FORMAT_VERSION,
                          "Bridge checkpoint format version")
            producer = manifest.get("bridge_version")
            if not isinstance(producer, str) or not producer.strip():
                raise ValueError("Checkpoint is missing a valid Bridge producer version.")
            if manifest.get("rebound_version") != rebound.__version__:
                raise ValueError(
                    f"Checkpoint REBOUND version {manifest.get('rebound_version')!r} "
                    f"does not match installed version {rebound.__version__!r}."
                )
            if manifest.get("abi") != _abi():
                raise ValueError("Checkpoint platform/ABI does not match this runtime.")
            metadata = manifest.get("metadata")
            entries = manifest.get("simulations")
            if not isinstance(metadata, dict) or not isinstance(entries, list) or not entries:
                raise ValueError("Checkpoint is missing simulation or hierarchy metadata.")
            if require_bridge_schema or "bridge_schema" in metadata:
                _check_schema(metadata.get("bridge_schema"), _BRIDGE_SCHEMA_VERSION,
                              "Bridge checkpoint schema")
            if "physical_tree" in metadata:
                tree = metadata["physical_tree"]
                if not isinstance(tree, dict):
                    raise ValueError("Invalid physical-tree schema metadata.")
                _check_schema(tree.get("schema"), _TREE_SCHEMA_VERSION,
                              "physical-tree schema")
            component_names = set()
            expected_members = {"manifest.json"}
            # 全部组件校验通过后再恢复，摘要不代表来源可信。
            for index, entry in enumerate(entries):
                if not isinstance(entry, dict):
                    raise ValueError("Invalid checkpoint simulation descriptor.")
                name = entry.get("name")
                member = entry.get("member")
                if not isinstance(name, str) or not name or name in component_names:
                    raise ValueError("Invalid or duplicate checkpoint component name.")
                if member != f"simulations/{index}.bin":
                    raise ValueError("Invalid checkpoint simulation member.")
                component_names.add(name)
                expected_members.add(member)
                info = archive.getinfo(member)
                if info.file_size != entry.get("size"):
                    raise ValueError(f"Checkpoint size mismatch for component {name!r}.")
                payload = archive.read(member)
                if hashlib.sha256(payload).hexdigest() != entry.get("sha256"):
                    raise ValueError(f"Checkpoint checksum mismatch for component {name!r}.")
                if not payload.startswith(b"REBOUND Binary File."):
                    raise ValueError(f"Invalid REBOUND archive for component {name!r}.")
                payloads.append((name, payload))
            if set(names) != expected_members:
                raise ValueError("Checkpoint has unexpected archive members.")
    except (zipfile.BadZipFile, KeyError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"Cannot read checkpoint archive: {exc}") from exc

    simulations = {}
    try:
        with warnings.catch_warnings():
            # 恢复警告按失败处理，避免带着缺失状态继续积分。
            warnings.simplefilter("error", RuntimeWarning)
            for name, payload in payloads:
                simulations[name] = rebound.Simulation(payload)
    except (RuntimeError, RuntimeWarning) as exc:
        raise ValueError(f"REBOUND could not restore checkpoint: {exc}") from exc
    return simulations, metadata
