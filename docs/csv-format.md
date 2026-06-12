# CSV Format Guide

KinetiKit uses **tidy CSVs**: each row is one observation, numeric values live in their own columns, and optional metadata columns can travel along without breaking analysis.

## Enzyme initial-rate CSV

Use this when you already have initial rates from any assay type.

Required columns, using any accepted name:

| Concept | Accepted names |
|---|---|
| substrate concentration | `substrate`, `substrate_uM`, `substrate_um`, `S`, `concentration` |
| initial rate | `rate`, `velocity`, `v0`, `vo`, `initial_rate`, `initial_rate_uM_per_s`, `initial_rate_mM_per_s` |

Optional columns:

- `replicate`
- `enzyme_uM` or another enzyme concentration metadata field
- `condition`
- `experiment_id`
- `temperature_C`
- `pH`
- `notes`

Duplicate substrate rows are averaged before fitting.

## General chemistry integrated rate-law CSV

Use this for concentration-vs-time data.

| Concept | Accepted names |
|---|---|
| time | `time`, `time_s`, `time_sec`, `seconds`, `t` |
| concentration | `concentration`, `conc`, `reactant`, `reactant_concentration`, `reactant_concentration_mM` |

KinetiKit can fit order `0`, `1`, `2`, or `auto`.

## Method-of-initial-rates CSV

Use this when you vary reactant concentrations and measure initial rates.

Required style:

- one or more columns beginning with `reactant_`
- one rate column such as `initial_rate_mM_per_s`

Example:

```csv
reactant_A_mM,reactant_B_mM,initial_rate_mM_per_s
10,5,0.55
20,5,1.10
10,10,2.20
20,10,4.40
```

KinetiKit fits:

```text
rate = k * A^m_A * B^m_B * ...
```

and reports `k`, `order_reactant_A_mM`, `order_reactant_B_mM`, etc.
