import math
from pathlib import Path

import pandas as pd
import pytest

from kinetikit import (
    fit_michaelis_menten,
    fit_substrate_inhibition,
    fit_rate_law,
    analyze_enzyme_csv,
    analyze_rate_csv,
)


def test_fit_michaelis_menten_recovers_parameters_and_derived_metrics():
    substrate = [2.5, 5, 10, 20, 40, 80, 160]
    vmax = 0.84
    km = 12.5
    rates = [vmax * s / (km + s) for s in substrate]

    result = fit_michaelis_menten(substrate, rates, enzyme_conc=0.25)

    assert result.model == "michaelis_menten"
    assert result.n_points == 7
    assert result.parameters["vmax"] == pytest.approx(vmax, rel=1e-3)
    assert result.parameters["km"] == pytest.approx(km, rel=1e-3)
    assert result.derived["kcat"] == pytest.approx(vmax / 0.25, rel=1e-3)
    assert result.derived["catalytic_efficiency"] == pytest.approx((vmax / 0.25) / km, rel=1e-3)
    assert result.goodness["r_squared"] > 0.9999
    assert "v = vmax * [S] / (km + [S])" in result.equation


def test_fit_substrate_inhibition_recovers_ki_for_high_substrate_dropoff():
    substrate = [2.5, 5, 10, 20, 40, 80, 160, 320]
    vmax = 1.2
    km = 9.0
    ki = 140.0
    rates = [vmax * s / (km + s * (1 + s / ki)) for s in substrate]

    result = fit_substrate_inhibition(substrate, rates, enzyme_conc=0.3)

    assert result.model == "substrate_inhibition"
    assert result.parameters["vmax"] == pytest.approx(vmax, rel=1e-2)
    assert result.parameters["km"] == pytest.approx(km, rel=2e-2)
    assert result.parameters["ki"] == pytest.approx(ki, rel=2e-2)
    assert result.derived["kcat"] == pytest.approx(vmax / 0.3, rel=1e-2)
    assert result.goodness["r_squared"] > 0.999


def test_fit_rate_law_identifies_second_order_from_concentration_time_course():
    k = 0.18
    a0 = 0.1
    time = [0, 5, 10, 15, 20, 25, 30]
    concentration = [1 / (1 / a0 + k * t) for t in time]

    result = fit_rate_law(time, concentration, order="auto")

    assert result.model == "second_order"
    assert result.parameters["k"] == pytest.approx(k, rel=1e-3)
    assert result.parameters["initial_concentration"] == pytest.approx(a0, rel=1e-3)
    assert result.goodness["r_squared"] > 0.9999
    assert result.equation == "1/[A] = kt + 1/[A]0"


def test_csv_analyzers_accept_generic_column_names(tmp_path):
    enzyme_csv = tmp_path / "enzyme_rates.csv"
    pd.DataFrame(
        {
            "substrate": [2.5, 5, 10, 20, 40, 80, 160],
            "rate": [0.84 * s / (12.5 + s) for s in [2.5, 5, 10, 20, 40, 80, 160]],
        }
    ).to_csv(enzyme_csv, index=False)

    enzyme_result = analyze_enzyme_csv(enzyme_csv, model="michaelis-menten", enzyme_conc=0.25)
    assert enzyme_result.parameters["km"] == pytest.approx(12.5, rel=1e-3)

    rate_csv = tmp_path / "rate_course.csv"
    k = 0.18
    a0 = 0.1
    times = [0, 5, 10, 15, 20, 25, 30]
    pd.DataFrame(
        {
            "time": times,
            "concentration": [1 / (1 / a0 + k * t) for t in times],
        }
    ).to_csv(rate_csv, index=False)

    rate_result = analyze_rate_csv(rate_csv, order="auto")
    assert rate_result.model == "second_order"
    assert rate_result.parameters["k"] == pytest.approx(k, rel=1e-3)
