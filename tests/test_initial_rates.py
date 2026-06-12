import pandas as pd
import pytest

from kinetikit import fit_initial_rate_orders, analyze_initial_rates_csv


def test_fit_initial_rate_orders_recovers_power_law_orders_for_two_reactants():
    rows = []
    for a in [10, 20, 40, 80]:
        for b in [5, 10, 20]:
            rows.append((a, b, 0.0022 * (a ** 1.0) * (b ** 2.0)))

    result = fit_initial_rate_orders(
        reactants={"A": [row[0] for row in rows], "B": [row[1] for row in rows]},
        rates=[row[2] for row in rows],
    )

    assert result.model == "power_law_initial_rates"
    assert result.parameters["order_A"] == pytest.approx(1.0, rel=1e-3)
    assert result.parameters["order_B"] == pytest.approx(2.0, rel=1e-3)
    assert result.parameters["k"] == pytest.approx(0.0022, rel=1e-3)
    assert result.goodness["r_squared"] > 0.9999
    assert result.equation == "rate = k * A^m_A * B^m_B"


def test_analyze_initial_rates_csv_detects_reactant_columns(tmp_path):
    csv_path = tmp_path / "two_reactants.csv"
    rows = []
    for a in [10, 20, 40, 80]:
        for b in [5, 10, 20]:
            rows.append({"reactant_A_mM": a, "reactant_B_mM": b, "initial_rate_mM_per_s": 0.0022 * a * b**2})
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    result = analyze_initial_rates_csv(csv_path)

    assert result.parameters["order_reactant_A_mM"] == pytest.approx(1.0, rel=1e-3)
    assert result.parameters["order_reactant_B_mM"] == pytest.approx(2.0, rel=1e-3)
    assert result.parameters["k"] == pytest.approx(0.0022, rel=1e-3)
