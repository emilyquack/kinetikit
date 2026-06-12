from __future__ import annotations

import numpy as np


def michaelis_menten(substrate, vmax: float, km: float):
    s = np.asarray(substrate, dtype=float)
    return vmax * s / (km + s)


def substrate_inhibition(substrate, vmax: float, km: float, ki: float):
    """Substrate inhibition model: v = Vmax*S / (Km + S*(1 + S/Ki))."""
    s = np.asarray(substrate, dtype=float)
    return vmax * s / (km + s * (1 + s / ki))


def zero_order(time, a0: float, k: float):
    t = np.asarray(time, dtype=float)
    return a0 - k * t


def first_order(time, a0: float, k: float):
    t = np.asarray(time, dtype=float)
    return a0 * np.exp(-k * t)


def second_order(time, a0: float, k: float):
    t = np.asarray(time, dtype=float)
    return 1 / (1 / a0 + k * t)
