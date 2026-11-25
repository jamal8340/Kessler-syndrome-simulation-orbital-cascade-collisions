"""Kessler simulation engine: analytic orbits, breakups, close approaches, scenarios.

The engine has no networking; the Unity server (step C) is built on top of
``kessler_sim.engine.Simulation``.
"""

from .breakup import BreakupConfig, BreakupEvent, BreakupParams
from .engine import Simulation, Snapshot
from .scenario import available_scenarios, load_scenario

__all__ = ["BreakupConfig", "BreakupEvent", "BreakupParams", "Simulation", "Snapshot",
           "available_scenarios", "load_scenario"]
