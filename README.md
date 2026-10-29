# State-level CCEMG versus ABC whole-sum DWB (`l_N = 1`)

This package estimates the Section 5 wheat-yield model separately for New South Wales (NSW), South Australia (SA), and Victoria (VIC). It compares standard CCEMG with the analytical bias-corrected estimator (ABC); the jackknife estimator is omitted.

The supplied balanced panel covers 1982–2016 (`T = 35`). NSW has 36 regions, SA has 25, and VIC has 24. The dependent variable is `ln_yields`; the regressors are `FD`, `DDlow`, `DDmod`, `DDhigh`, `RF`, `FD_RF`, and `DDhigh_RF`.

## Estimation and inference

Each individual CCE regression includes an intercept, the seven regressors, and the state-specific cross-sectional averages of log yield and the seven regressors. The CCEMG estimate is the mean of the individual CCE slopes. Its standard covariance estimator is

```text
V_CCEMG = [N(N-1)]^(-1) sum_i (beta_i - beta_MG)(beta_i - beta_MG)'.
```

The ABC estimator uses the corrected subtraction sign from Equation (44):

```text
beta_ABC = beta_MG - (b0_hat + b1_hat)/sqrt(NT).
```

The corresponding corrected individual slopes are

```text
beta_ABC,i = beta_i - (b0_i + b1_i)/N,
```

whose mean is `beta_ABC`. The CCEMG regressions retain all cross-sectional averages. The nuisance factor space used only for the analytical correction defaults to the first two principal components (`FACTOR_COUNT=2`).

For ABC inference, the code implements the covariance formula with `l_N = 1`. Let

```text
Z_i   = sqrt(T) (beta_ABC,i - beta_ABC),
Z_i^c = sqrt[N/(N-1)] Z_i.
```

Since `l_N = 1`, the kernel removes cross-unit terms and the conditional DWB covariance is

```text
Sigma_DWB = N^(-1) sum_i Z_i^c Z_i^c'.
```

Thus the covariance of `beta_ABC` used for its standard errors is exactly

```text
V_ABC = Sigma_DWB/(NT)
      = [N(N-1)]^(-1) sum_i
        (beta_ABC,i - beta_ABC)(beta_ABC,i - beta_ABC)'.
```

Bootstrap statistics use independent cross-sectional standard-normal multipliers. The asymmetric basic interval is

```text
[ beta_ABC - q*_(1-alpha/2)/sqrt(NT),
  beta_ABC - q*_(alpha/2)/sqrt(NT) ].
```

The default run uses 9,999 bootstrap replications, 95% confidence, and seed `20261001`.

## Run

Python:

```bash
bash run_empirical_analysis_python.sh
```

R (requires `MASS`, `readxl`, and `ggplot2`):

```bash
bash run_empirical_analysis.sh
```

Configuration variables are `BOOT_REPS`, `CONFIDENCE_LEVEL`, `FACTOR_COUNT`, `EMPIRICAL_SEED`, `DATA_FILE`, and `OUTPUT_DIR`. The Python implementation exposes matching command-line options. The historical script filenames retain `hybrid` only for drop-in compatibility; their implemented ABC inference is the whole-sum DWB with `l_N=1` described above.

## Outputs

- `empirical_ccemg_abc_results.csv`: estimates, standard errors, interval endpoints, interval lengths, and normal-reference p-values.
- `confidence_interval_lengths.csv`: direct CCEMG/ABC interval-length comparison.
- `individual_cce_slopes_and_abc_components.csv`: unit CCE slopes, correction components, corrected ABC slopes, and centered DWB contributions.
- `empirical_diagnostics.csv`: panel/rank diagnostics and numerical verification of the `l_N=1` covariance identity.
- `NSW_mean_effects_with_shaded_CI.png`, `SA_mean_effects_with_shaded_CI.png`, and `VIC_mean_effects_with_shaded_CI.png`.
- `RESULTS_SUMMARY.md` and `run_configuration.json`.
