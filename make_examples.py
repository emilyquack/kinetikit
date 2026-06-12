from pathlib import Path
import csv
import math
import random
import json
import zipfile

base = Path(__file__).resolve().parent
examples = base / 'examples'
examples.mkdir(parents=True, exist_ok=True)
random.seed(42)

# 1) Generic enzyme Michaelis-Menten initial rates (replicate rows, no plate-reader assumptions)
S_values = [1, 2.5, 5, 10, 20, 40, 80, 160, 320]
Vmax, Km = 1.25, 18.0  # µM/s, µM
rows = []
for S in S_values:
    true_v = Vmax * S / (Km + S)
    for rep in range(1, 4):
        noise = random.gauss(0, true_v * 0.035 + 0.004)
        rows.append({
            'experiment_id': 'MM_demo_001',
            'condition': 'control',
            'replicate': rep,
            'substrate_uM': S,
            'initial_rate_uM_per_s': round(max(true_v + noise, 0), 5),
            'enzyme_uM': 0.15,
            'temperature_C': 25,
            'pH': 7.5,
            'notes': 'Generic purified enzyme initial-rate example; rates can come from any assay type.'
        })
with (examples / 'enzyme_michaelis_menten_initial_rates.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# 2) Enzyme substrate inhibition initial rates (higher [S] starts reducing velocity)
S_values2 = [1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1000]
Vmax, Km, Ki = 1.55, 16.0, 420.0
rows = []
for S in S_values2:
    true_v = Vmax * S / (Km + S + (S * S / Ki))
    for rep in range(1, 4):
        noise = random.gauss(0, true_v * 0.04 + 0.004)
        rows.append({
            'experiment_id': 'SUBINHIB_demo_001',
            'condition': 'control',
            'replicate': rep,
            'substrate_uM': S,
            'initial_rate_uM_per_s': round(max(true_v + noise, 0), 5),
            'enzyme_uM': 0.15,
            'temperature_C': 25,
            'pH': 7.5,
            'expected_model': 'substrate_inhibition',
            'notes': 'Example where very high substrate inhibits the reaction; useful for comparing MM vs substrate-inhibition fits.'
        })
with (examples / 'enzyme_substrate_inhibition_initial_rates.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# 3) Enzyme inhibition screen: rates at several inhibitor concentrations
S_values3 = [5, 10, 20, 40, 80, 160]
inhibitors = [0, 5, 20, 80]
Vmax, Km, Ki_comp = 1.20, 14.0, 22.0
rows = []
for I in inhibitors:
    for S in S_values3:
        true_v = Vmax * S / (Km * (1 + I / Ki_comp) + S)
        for rep in range(1, 3):
            noise = random.gauss(0, true_v * 0.04 + 0.003)
            rows.append({
                'experiment_id': 'INHIB_demo_001',
                'condition': f'inhibitor_{I}_uM',
                'replicate': rep,
                'substrate_uM': S,
                'inhibitor_uM': I,
                'initial_rate_uM_per_s': round(max(true_v + noise, 0), 5),
                'enzyme_uM': 0.12,
                'temperature_C': 25,
                'pH': 7.5,
                'expected_model': 'competitive_inhibition',
                'notes': 'Generic inhibition example with multiple inhibitor concentrations.'
            })
with (examples / 'enzyme_inhibition_initial_rates.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# 4) Gen chem integrated rate law concentration-vs-time example (non-enzyme)
rows = []
A0, k = 100.0, 0.0045  # mM, s^-1
for t in range(0, 1001, 100):
    true_A = A0 * math.exp(-k * t)
    for rep in range(1, 4):
        noise = random.gauss(0, 0.7)
        rows.append({
            'reaction_id': 'GENCHEM_first_order_demo',
            'condition': '25C_buffered',
            'replicate': rep,
            'time_s': t,
            'reactant_concentration_mM': round(max(true_A + noise, 0), 4),
            'temperature_C': 25,
            'expected_order': 1,
            'notes': 'Non-enzyme general chemistry concentration-vs-time data for integrated rate-law fitting.'
        })
with (examples / 'genchem_first_order_timecourse.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# 5) Gen chem initial-rate method example (vary reactant concentrations, record initial rate)
rows = []
for A in [10, 20, 40, 80]:
    for B in [5, 10, 20]:
        true_rate = 0.0022 * (A ** 1) * (B ** 2)  # designed first-order in A, second-order in B
        for rep in range(1, 4):
            noise = random.gauss(0, true_rate * 0.035)
            rows.append({
                'reaction_id': 'GENCHEM_initial_rates_demo',
                'condition': 'method_of_initial_rates',
                'replicate': rep,
                'reactant_A_mM': A,
                'reactant_B_mM': B,
                'initial_rate_mM_per_s': round(max(true_rate + noise, 0), 5),
                'temperature_C': 25,
                'notes': 'Designed so rate is approximately first order in A and second order in B.'
            })
with (examples / 'genchem_initial_rates_two_reactants.csv').open('w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

readme = '''# Kinetikit example CSVs

These are starter datasets for **Kinetikit**, a generic kinetics analysis toolkit. They intentionally avoid plate-reader-specific raw absorbance formats so the package can work with rates or concentration-time data from many assay types.

## Files

### `enzyme_michaelis_menten_initial_rates.csv`
Generic enzyme initial-rate data for a Michaelis-Menten fit.

Required-ish columns for this style:
- `substrate_uM` — substrate concentration
- `initial_rate_uM_per_s` — initial velocity
- `replicate` — replicate number

Helpful optional columns:
- `enzyme_uM` — enzyme concentration, useful for kcat = Vmax / [E]
- `condition`, `experiment_id`, `temperature_C`, `pH`, `notes`

### `enzyme_substrate_inhibition_initial_rates.csv`
Enzyme initial-rate data where high substrate concentrations reduce rate. This is for fitting/comparing:

`v = Vmax * S / (Km + S + S^2/Ki)`

It has the same core columns as the Michaelis-Menten example plus `expected_model`.

### `enzyme_inhibition_initial_rates.csv`
Initial-rate data measured at several inhibitor concentrations. Useful for competitive/noncompetitive/uncompetitive inhibition workflows.

Extra column:
- `inhibitor_uM` — inhibitor concentration

### `genchem_first_order_timecourse.csv`
Non-enzyme general chemistry concentration-vs-time data. Useful for integrated rate-law analysis.

Core columns:
- `time_s`
- `reactant_concentration_mM`
- `replicate`

### `genchem_initial_rates_two_reactants.csv`
General chemistry method-of-initial-rates dataset with two reactants.

Core columns:
- `reactant_A_mM`
- `reactant_B_mM`
- `initial_rate_mM_per_s`
- `replicate`

## Design principle

Kinetikit should accept tidy CSVs: one observation per row, numeric columns for concentrations/rates/time, and optional metadata columns for condition labels.
'''
(examples / 'README.md').write_text(readme)

zip_path = base / 'kinetikit_example_csvs.zip'
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for path in sorted(examples.iterdir()):
        z.write(path, arcname=f'examples/{path.name}')

print(json.dumps({
    'base': str(base),
    'examples': [p.name for p in sorted(examples.glob('*.csv'))],
    'readme': str(examples / 'README.md'),
    'zip': str(zip_path),
}, indent=2))
