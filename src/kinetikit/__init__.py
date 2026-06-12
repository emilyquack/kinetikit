"""KinetiKit: friendly enzyme kinetics and chemistry rate-law fitting."""

from .fitting import (
    analyze_enzyme_csv,
    analyze_initial_rates_csv,
    analyze_rate_csv,
    fit_initial_rate_orders,
    fit_michaelis_menten,
    fit_rate_law,
    fit_substrate_inhibition,
)
from .results import FitResult

__all__ = [
    "FitResult",
    "fit_michaelis_menten",
    "fit_substrate_inhibition",
    "fit_rate_law",
    "fit_initial_rate_orders",
    "analyze_enzyme_csv",
    "analyze_rate_csv",
    "analyze_initial_rates_csv",
]

__version__ = "0.1.0"
