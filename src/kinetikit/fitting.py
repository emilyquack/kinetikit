from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .models import michaelis_menten, substrate_inhibition
from .results import FitResult

SUBSTRATE_ALIASES = ("substrate", "substrate_um", "substrate_uM", "s", "S", "concentration", "concentration_um")
RATE_ALIASES = (
    "rate",
    "velocity",
    "v0",
    "vo",
    "initial_rate",
    "initial_velocity",
    "initial_rate_um_per_s",
    "initial_rate_uM_per_s",
    "initial_rate_mm_per_s",
    "initial_rate_mM_per_s",
)
TIME_ALIASES = ("time", "time_s", "time_sec", "seconds", "t")
CONC_ALIASES = (
    "concentration",
    "conc",
    "a",
    "A",
    "reactant",
    "reactant_concentration",
    "reactant_concentration_mm",
    "reactant_concentration_mM",
)


def _as_float_array(values: Iterable[float], name: str) -> np.ndarray:
    arr = np.asarray(list(values), dtype=float)
    if arr.ndim != 1 or arr.size == 0:
        raise ValueError(f"{name} must be a non-empty one-dimensional sequence")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values")
    return arr


def _r_squared(observed: np.ndarray, predicted: np.ndarray) -> float:
    ss_res = float(np.sum((observed - predicted) ** 2))
    ss_tot = float(np.sum((observed - np.mean(observed)) ** 2))
    if ss_tot == 0:
        return 1.0 if ss_res == 0 else float("nan")
    return 1 - ss_res / ss_tot


def _rmse(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.sqrt(np.mean((observed - predicted) ** 2)))


def _clean_xy(x, y, x_name: str, y_name: str) -> tuple[np.ndarray, np.ndarray]:
    x_arr = _as_float_array(x, x_name)
    y_arr = _as_float_array(y, y_name)
    if x_arr.size != y_arr.size:
        raise ValueError(f"{x_name} and {y_name} must have the same length")
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    x_arr, y_arr = x_arr[mask], y_arr[mask]
    if x_arr.size < 3:
        raise ValueError("At least three valid points are required for fitting")
    if np.any(x_arr < 0):
        raise ValueError(f"{x_name} must be non-negative")
    return x_arr, y_arr


def _parameter_errors(covariance: np.ndarray, names: list[str]) -> dict[str, float]:
    if covariance.size == 0 or not np.all(np.isfinite(covariance)):
        return {f"{name}_stderr": float("nan") for name in names}
    return {f"{name}_stderr": float(np.sqrt(max(covariance[i, i], 0))) for i, name in enumerate(names)}


def fit_michaelis_menten(substrate, rates, enzyme_conc: float | None = None) -> FitResult:
    s, v = _clean_xy(substrate, rates, "substrate", "rates")
    if np.any(s <= 0):
        raise ValueError("Michaelis-Menten fitting requires positive substrate concentrations")
    if np.any(v < 0):
        raise ValueError("Rates must be non-negative")

    vmax_guess = float(max(v))
    half_max = vmax_guess / 2 if vmax_guess else 1.0
    km_guess = float(s[np.argmin(np.abs(v - half_max))]) if s.size else 1.0
    bounds = ([0, 0], [np.inf, np.inf])
    params, cov = curve_fit(michaelis_menten, s, v, p0=[vmax_guess, km_guess], bounds=bounds, maxfev=10000)
    vmax, km = [float(x) for x in params]
    pred = michaelis_menten(s, vmax, km)
    derived: dict[str, float] = {}
    warnings: list[str] = []
    if enzyme_conc is not None:
        if enzyme_conc <= 0:
            raise ValueError("enzyme_conc must be positive when provided")
        derived["kcat"] = vmax / enzyme_conc
        derived["catalytic_efficiency"] = derived["kcat"] / km if km else float("nan")
    else:
        warnings.append("enzyme_conc was not provided, so kcat and catalytic efficiency were not calculated")

    parameter_errors = _parameter_errors(cov, ["vmax", "km"])
    return FitResult(
        model="michaelis_menten",
        equation="v = vmax * [S] / (km + [S])",
        parameters={"vmax": vmax, "km": km, **parameter_errors},
        derived=derived,
        goodness={"r_squared": _r_squared(v, pred), "rmse": _rmse(v, pred)},
        n_points=int(s.size),
        warnings=warnings,
        data={"substrate": s.tolist(), "rate": v.tolist()},
        predictions=[float(x) for x in pred],
    )


def fit_substrate_inhibition(substrate, rates, enzyme_conc: float | None = None) -> FitResult:
    s, v = _clean_xy(substrate, rates, "substrate", "rates")
    if np.any(s <= 0):
        raise ValueError("Substrate inhibition fitting requires positive substrate concentrations")
    if np.any(v < 0):
        raise ValueError("Rates must be non-negative")

    vmax_guess = float(max(v) * 1.2)
    km_guess = float(np.median(s))
    ki_guess = float(max(s))
    params, cov = curve_fit(
        substrate_inhibition,
        s,
        v,
        p0=[vmax_guess, km_guess, ki_guess],
        bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
        maxfev=20000,
    )
    vmax, km, ki = [float(x) for x in params]
    pred = substrate_inhibition(s, vmax, km, ki)
    derived: dict[str, float] = {}
    warnings: list[str] = []
    if enzyme_conc is not None:
        if enzyme_conc <= 0:
            raise ValueError("enzyme_conc must be positive when provided")
        derived["kcat"] = vmax / enzyme_conc
        derived["catalytic_efficiency"] = derived["kcat"] / km if km else float("nan")
    else:
        warnings.append("enzyme_conc was not provided, so kcat and catalytic efficiency were not calculated")

    parameter_errors = _parameter_errors(cov, ["vmax", "km", "ki"])
    return FitResult(
        model="substrate_inhibition",
        equation="v = vmax * [S] / (km + [S] * (1 + [S]/ki))",
        parameters={"vmax": vmax, "km": km, "ki": ki, **parameter_errors},
        derived=derived,
        goodness={"r_squared": _r_squared(v, pred), "rmse": _rmse(v, pred)},
        n_points=int(s.size),
        warnings=warnings,
        data={"substrate": s.tolist(), "rate": v.tolist()},
        predictions=[float(x) for x in pred],
    )


def _linear_fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float, np.ndarray, float]:
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    return float(slope), float(intercept), pred, _r_squared(y, pred)


def _fit_rate_order(time: np.ndarray, concentration: np.ndarray, order: int) -> FitResult:
    if order == 0:
        slope, intercept, transformed_pred, r2 = _linear_fit(time, concentration)
        k = -slope
        pred = intercept - k * time
        return FitResult(
            model="zero_order",
            equation="[A] = [A]0 - kt",
            parameters={"k": float(k), "initial_concentration": float(intercept)},
            goodness={"r_squared": r2, "rmse": _rmse(concentration, pred)},
            n_points=int(time.size),
            data={"time": time.tolist(), "concentration": concentration.tolist()},
            predictions=[float(x) for x in pred],
        )
    if order == 1:
        if np.any(concentration <= 0):
            raise ValueError("First-order fitting requires positive concentrations")
        slope, intercept, transformed_pred, r2 = _linear_fit(time, np.log(concentration))
        k = -slope
        a0 = float(np.exp(intercept))
        pred = a0 * np.exp(-k * time)
        return FitResult(
            model="first_order",
            equation="ln[A] = -kt + ln[A]0",
            parameters={"k": float(k), "initial_concentration": a0},
            goodness={"r_squared": r2, "rmse": _rmse(concentration, pred)},
            n_points=int(time.size),
            data={"time": time.tolist(), "concentration": concentration.tolist()},
            predictions=[float(x) for x in pred],
        )
    if order == 2:
        if np.any(concentration <= 0):
            raise ValueError("Second-order fitting requires positive concentrations")
        slope, intercept, transformed_pred, r2 = _linear_fit(time, 1 / concentration)
        k = slope
        a0 = 1 / intercept
        pred = 1 / (1 / a0 + k * time)
        return FitResult(
            model="second_order",
            equation="1/[A] = kt + 1/[A]0",
            parameters={"k": float(k), "initial_concentration": float(a0)},
            goodness={"r_squared": r2, "rmse": _rmse(concentration, pred)},
            n_points=int(time.size),
            data={"time": time.tolist(), "concentration": concentration.tolist()},
            predictions=[float(x) for x in pred],
        )
    raise ValueError("order must be 0, 1, 2, or 'auto'")


def fit_initial_rate_orders(reactants: dict[str, Iterable[float]], rates: Iterable[float]) -> FitResult:
    """Fit method-of-initial-rates data to rate = k * A^m * B^n * ... ."""
    if not reactants:
        raise ValueError("At least one reactant concentration series is required")
    names = list(reactants.keys())
    arrays = [_as_float_array(reactants[name], name) for name in names]
    rate_arr = _as_float_array(rates, "rates")
    if any(arr.size != rate_arr.size for arr in arrays):
        raise ValueError("All reactant arrays and rates must have the same length")
    if any(np.any(arr <= 0) for arr in arrays):
        raise ValueError("Initial-rate order fitting requires positive reactant concentrations")
    if np.any(rate_arr <= 0):
        raise ValueError("Initial-rate order fitting requires positive rates")

    design = np.column_stack([np.ones(rate_arr.size), *[np.log(arr) for arr in arrays]])
    coefficients, *_ = np.linalg.lstsq(design, np.log(rate_arr), rcond=None)
    ln_k = float(coefficients[0])
    orders = [float(x) for x in coefficients[1:]]
    k = float(np.exp(ln_k))
    pred_log = design @ coefficients
    pred = np.exp(pred_log)
    parameters = {"k": k}
    parameters.update({f"order_{name}": order for name, order in zip(names, orders)})
    symbolic_terms = [f"{name}^m_{name}" for name in names]
    return FitResult(
        model="power_law_initial_rates",
        equation=f"rate = k * {' * '.join(symbolic_terms)}",
        parameters=parameters,
        goodness={"r_squared": _r_squared(np.log(rate_arr), pred_log), "rmse": _rmse(rate_arr, pred)},
        n_points=int(rate_arr.size),
        data={**{name: arr.tolist() for name, arr in zip(names, arrays)}, "rate": rate_arr.tolist()},
        predictions=[float(x) for x in pred],
    )


def fit_rate_law(time, concentration, order: int | str = "auto") -> FitResult:
    t, c = _clean_xy(time, concentration, "time", "concentration")
    if np.any(c < 0):
        raise ValueError("concentration must be non-negative")
    if order == "auto":
        candidates = [_fit_rate_order(t, c, candidate) for candidate in (0, 1, 2)]
        best = max(candidates, key=lambda result: result.goodness["r_squared"])
        best.derived["tested_orders"] = 3
        return best
    return _fit_rate_order(t, c, int(order))


def _find_column(df: pd.DataFrame, aliases: tuple[str, ...], role: str) -> str:
    normalized = {str(col).strip().lower(): col for col in df.columns}
    for alias in aliases:
        key = alias.strip().lower()
        if key in normalized:
            return normalized[key]
    raise ValueError(f"Could not find a {role} column. Tried: {', '.join(aliases)}")


def _read_csv(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def analyze_enzyme_csv(path: str | Path, model: str = "michaelis-menten", enzyme_conc: float | None = None) -> FitResult:
    df = _read_csv(path)
    substrate_col = _find_column(df, SUBSTRATE_ALIASES, "substrate")
    rate_col = _find_column(df, RATE_ALIASES, "rate")
    clean = df[[substrate_col, rate_col]].dropna().copy()
    # Average duplicate substrate concentrations, a common kinetics workflow.
    grouped = clean.groupby(substrate_col, as_index=False)[rate_col].mean().sort_values(substrate_col)
    model_key = model.lower().replace("_", "-")
    if model_key in {"michaelis-menten", "mm"}:
        return fit_michaelis_menten(grouped[substrate_col], grouped[rate_col], enzyme_conc=enzyme_conc)
    if model_key in {"substrate-inhibition", "substrate inhibition", "inhibition"}:
        return fit_substrate_inhibition(grouped[substrate_col], grouped[rate_col], enzyme_conc=enzyme_conc)
    raise ValueError("model must be 'michaelis-menten' or 'substrate-inhibition'")


def analyze_initial_rates_csv(path: str | Path, reactant_columns: list[str] | None = None) -> FitResult:
    df = _read_csv(path)
    rate_col = _find_column(df, RATE_ALIASES, "rate")
    if reactant_columns is None:
        reactant_columns = [
            col for col in df.columns
            if str(col).strip().lower().startswith("reactant_")
            and "concentration" not in str(col).strip().lower()
            and col != rate_col
        ]
    if not reactant_columns:
        raise ValueError("Could not detect reactant columns. Use columns named like reactant_A_mM and reactant_B_mM.")
    missing = [col for col in reactant_columns if col not in df.columns]
    if missing:
        raise ValueError(f"Reactant columns not found: {', '.join(missing)}")
    clean = df[[*reactant_columns, rate_col]].dropna().copy()
    grouped = clean.groupby(reactant_columns, as_index=False)[rate_col].mean()
    return fit_initial_rate_orders(
        reactants={col: grouped[col] for col in reactant_columns},
        rates=grouped[rate_col],
    )


def analyze_rate_csv(path: str | Path, order: int | str = "auto") -> FitResult:
    df = _read_csv(path)
    time_col = _find_column(df, TIME_ALIASES, "time")
    conc_col = _find_column(df, CONC_ALIASES, "concentration")
    clean = df[[time_col, conc_col]].dropna().sort_values(time_col)
    return fit_rate_law(clean[time_col], clean[conc_col], order=order)
