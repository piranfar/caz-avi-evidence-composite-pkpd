"""Unit-explicit continuous-infusion calculations.

No study-specific parameter values are embedded in this module.
"""

from __future__ import annotations

import numpy as np


def steady_state_concentration(
    infusion_rate_mg_per_h: np.ndarray | float,
    clearance_l_per_h: np.ndarray | float,
) -> np.ndarray:
    """Return total steady-state concentration in mg/L.

    Css = infusion rate / clearance. Inputs must use mg/h and L/h.
    """
    rate = np.asarray(infusion_rate_mg_per_h, dtype=float)
    clearance = np.asarray(clearance_l_per_h, dtype=float)
    if np.any(rate < 0):
        raise ValueError("Infusion rate must be non-negative.")
    if np.any(clearance <= 0):
        raise ValueError("Clearance must be strictly positive.")
    return rate / clearance


def unbound_concentration(
    total_concentration_mg_per_l: np.ndarray | float,
    free_fraction: np.ndarray | float,
) -> np.ndarray:
    """Convert total concentration to unbound concentration."""
    total = np.asarray(total_concentration_mg_per_l, dtype=float)
    fraction = np.asarray(free_fraction, dtype=float)
    if np.any(total < 0):
        raise ValueError("Concentration must be non-negative.")
    if np.any((fraction < 0) | (fraction > 1)):
        raise ValueError("Free fraction must be between 0 and 1.")
    return total * fraction
