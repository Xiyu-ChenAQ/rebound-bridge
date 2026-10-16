"""Multirate integration of fixed REBOUND hierarchies."""

from importlib.metadata import version

from ._native import BridgeError
from .bridge import Bridge
from .checkpoint import CheckpointError
from .hierarchy import Subsystem
from .model import Body, Simulation

__version__ = version("rebound-bridge")

__all__ = ["Body", "Simulation", "Bridge", "BridgeError", "CheckpointError", "Subsystem", "__version__"]
