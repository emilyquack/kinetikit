# Kinetikit example CSVs

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
