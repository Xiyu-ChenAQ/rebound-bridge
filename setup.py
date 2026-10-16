"""Build the Bridge layer; use the loaded Python REBOUND runtime at execution."""
import sys
from pathlib import Path

import rebound

from setuptools import Extension, setup

if rebound.__version__ != "4.6.0":
    raise RuntimeError("Build REBOUND_BRIDGE with the official rebound==4.6.0 dependency.")

symbols = (
    "reb_simulation_integrate",
    "reb_simulation_synchronize",
    "reb_simulation_free",
    "reb_particles_transform_jacobi_to_inertial_posvel",
)
macros = [(name, "rbpy_" + name) for name in symbols]
macros += [("REB_BRIDGE_BUILD_SHARED", "1"), ("REB_BRIDGE_BUILDING_LIBRARY", "1")]
macros.append(("REB_BRIDGE_PYTHON", "1"))
rebound_include = Path(rebound.__file__).parent
if "REB_INTEGRATOR_CUSTOM" not in (rebound_include / "rebound.h").read_text():
    macros.append(("REB_BRIDGE_REBOUND_NO_CUSTOM", "1"))

setup(ext_modules=[Extension(
    "rebound_bridge._bridge",
    sources=[
        "src/bridge.c", "src/bridge_advance.c", "src/bridge_eos.c",
        "src/bridge_kick.c", "src/bridge_whfast.c", "src/bridge_python.c",
    ],
    include_dirs=["src/include", str(rebound_include)],
    define_macros=macros,
    extra_compile_args=["/std:c11", "/O2", "/utf-8"] if sys.platform == "win32" else ["-std=c99", "-O3"],
    libraries=[] if sys.platform == "win32" else ["m"],
)])
