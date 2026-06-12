import json
from pathlib import Path

import pandas as pd

from kinetikit.cli import main


def test_cli_enzyme_writes_json_and_plot(tmp_path):
    csv_path = tmp_path / "enzyme.csv"
    pd.DataFrame({"substrate": [2.5, 5, 10, 20, 40, 80, 160], "rate": [0.84 * s / (12.5 + s) for s in [2.5, 5, 10, 20, 40, 80, 160]]}).to_csv(csv_path, index=False)
    json_path = tmp_path / "result.json"
    plot_path = tmp_path / "plot.png"

    exit_code = main(["enzyme", str(csv_path), "--enzyme-conc", "0.25", "--json", str(json_path), "--plot", str(plot_path)])

    assert exit_code == 0
    data = json.loads(json_path.read_text())
    assert data["model"] == "michaelis_menten"
    assert data["derived"]["kcat"] > 0
    assert plot_path.exists()
    assert plot_path.stat().st_size > 1000


def test_cli_rate_auto_writes_second_order_json(tmp_path):
    csv_path = tmp_path / "rate.csv"
    k = 0.18
    a0 = 0.1
    times = [0, 5, 10, 15, 20, 25, 30]
    pd.DataFrame({"time": times, "concentration": [1 / (1 / a0 + k * t) for t in times]}).to_csv(csv_path, index=False)
    json_path = tmp_path / "rate.json"

    exit_code = main(["rate", str(csv_path), "--order", "auto", "--json", str(json_path)])

    assert exit_code == 0
    data = json.loads(json_path.read_text())
    assert data["model"] == "second_order"
    assert data["parameters"]["k"] > 0


def test_cli_initial_rates_writes_power_law_json(tmp_path):
    csv_path = tmp_path / "initial.csv"
    rows = []
    for a in [10, 20, 40, 80]:
        for b in [5, 10, 20]:
            rows.append({"reactant_A_mM": a, "reactant_B_mM": b, "initial_rate_mM_per_s": 0.0022 * a * b**2})
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    json_path = tmp_path / "initial.json"

    exit_code = main(["initial-rates", str(csv_path), "--json", str(json_path)])

    assert exit_code == 0
    data = json.loads(json_path.read_text())
    assert data["model"] == "power_law_initial_rates"
    assert data["parameters"]["order_reactant_A_mM"] > 0.9
    assert data["parameters"]["order_reactant_B_mM"] > 1.9
