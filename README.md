# KinetiKit

**KinetiKit** is a friendly Python toolkit for fitting enzyme kinetics and general chemistry rate-law data from tidy CSV files. It is intentionally **not tied to plate-reader/NADPH/NADH raw formats**: users provide substrate/rate or time/concentration data from any assay source.

Built as a resume-worthy scientific software project: tested Python package, command-line interface, example datasets, plots, JSON outputs, and clear CSV schemas.

## What it can analyze

### Enzyme kinetics

- Michaelis-Menten initial-rate fits
- Substrate inhibition fits for high-substrate dropoff
- Optional `kcat` and catalytic efficiency when enzyme concentration is provided
- Duplicate substrate concentrations are averaged automatically

### General chemistry kinetics

- Integrated rate-law fitting for concentration-vs-time data
  - zero order
  - first order
  - second order
  - automatic order selection by linearized R²
- Method of initial rates for multiple reactants
  - fits `rate = k * A^m * B^n * ...`
  - estimates reaction orders and rate constant from tidy CSVs

## Install locally

From this folder:

```bash
python -m pip install -e '.[dev]'
```

Run tests:

```bash
python -m pytest -q
```

## CSV formats

### 1. Enzyme initial rates

Minimum columns:

```csv
substrate,rate
2.5,0.140
5,0.240
10,0.373
```

Also accepted:

- `substrate_uM`, `substrate_um`, `S`
- `initial_rate_uM_per_s`, `initial_rate_mM_per_s`, `velocity`, `v0`, `vo`
- optional metadata like `replicate`, `condition`, `temperature_C`, `pH`, `enzyme_uM`

Example files:

- `examples/enzyme_michaelis_menten.csv`
- `examples/enzyme_substrate_inhibition.csv`
- `examples/enzyme_michaelis_menten_initial_rates.csv`
- `examples/enzyme_substrate_inhibition_initial_rates.csv`

### 2. General chemistry time course

Minimum columns:

```csv
time,concentration
0,0.100
5,0.0526
10,0.0357
```

Also accepted:

- `time_s`, `seconds`, `t`
- `reactant_concentration_mM`, `conc`, `reactant`

Example files:

- `examples/rate_zero_order.csv`
- `examples/rate_first_order.csv`
- `examples/rate_second_order.csv`
- `examples/genchem_first_order_timecourse.csv`

### 3. Method of initial rates

Minimum style:

```csv
reactant_A_mM,reactant_B_mM,initial_rate_mM_per_s
10,5,0.55
20,5,1.10
10,10,2.20
```

Any columns beginning with `reactant_` are treated as reactants. The rate column can be `rate`, `initial_rate`, or `initial_rate_mM_per_s`.

Example file:

- `examples/genchem_initial_rates_two_reactants.csv`

## Command line examples

Fit a Michaelis-Menten dataset and make a plot/JSON output:

```bash
kinetikit enzyme examples/enzyme_michaelis_menten.csv \
  --model michaelis-menten \
  --enzyme-conc 0.25 \
  --json outputs/mm_result.json \
  --plot outputs/mm_fit.png
```

Fit substrate inhibition:

```bash
kinetikit enzyme examples/enzyme_substrate_inhibition.csv \
  --model substrate-inhibition \
  --enzyme-conc 0.3
```

Auto-detect zero/first/second order from a time course:

```bash
kinetikit rate examples/rate_second_order.csv --order auto --json outputs/rate_result.json
```

Estimate reaction orders from method-of-initial-rates data:

```bash
kinetikit initial-rates examples/genchem_initial_rates_two_reactants.csv \
  --json outputs/initial_rates_result.json
```

## Python API examples

```python
from kinetikit import fit_michaelis_menten, fit_rate_law

result = fit_michaelis_menten(
    substrate=[2.5, 5, 10, 20, 40, 80, 160],
    rates=[0.14, 0.24, 0.37, 0.52, 0.64, 0.72, 0.78],
    enzyme_conc=0.25,
)
print(result.summary())
```

```python
from kinetikit import fit_initial_rate_orders

result = fit_initial_rate_orders(
    reactants={"A": [10, 20, 10, 20], "B": [5, 5, 10, 10]},
    rates=[0.55, 1.10, 2.20, 4.40],
)
print(result.parameters)
```

## Project structure

```text
kinetikit/
├── src/kinetikit/
│   ├── fitting.py      # fitting algorithms and CSV analyzers
│   ├── models.py       # kinetics model equations
│   ├── plotting.py     # PNG plot exporters
│   ├── results.py      # serializable FitResult object
│   └── cli.py          # command-line interface
├── tests/              # pytest suite
├── examples/           # example CSV inputs
└── pyproject.toml      # package metadata
```

## Resume bullet idea

> Built **KinetiKit**, a tested Python package and CLI for enzyme kinetics and general chemistry rate-law analysis, using agent-assisted development with test-driven workflows; supports Michaelis-Menten, substrate inhibition, integrated rate laws, method-of-initial-rates order estimation, CSV ingestion, JSON export, and publication-ready plots.

## Roadmap ideas

- Confidence intervals for all models
- Model comparison/AIC for Michaelis-Menten vs substrate inhibition
- Inhibition models: competitive, noncompetitive, uncompetitive
- A Streamlit demo app or Jupyter tutorial notebook
- PyPI packaging/release workflow
