#!/usr/bin/env python3
"""Empirical CCEMG versus ABC whole-sum DWB analysis with l_N = 1.

The model follows Section 5 of the localized-parameter working paper.  Each
state is analysed as a separate balanced panel.  Standard CCEMG inference uses
Pesaran's cross-sectional dispersion estimator. ABC inference uses the complete
centered contribution of each analytically corrected individual slope. With
l_N = 1, the cross-sectional multipliers are independent and the exact
finite-sample-corrected conditional covariance is used for standard errors.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


VARIABLES = ["FD", "DDlow", "DDmod", "DDhigh", "RF", "FD_RF", "DDhigh_RF"]
DISPLAY_NAMES = {
    "FD": "FD",
    "DDlow": "DDlow",
    "DDmod": "DDmod",
    "DDhigh": "DDhigh",
    "RF": "RF",
    "FD_RF": "FD × RF",
    "DDhigh_RF": "DDhigh × RF",
}
STATES = ("NSW", "SA", "VIC")
ABC_METHOD = "ABC + whole-sum DWB"
METHOD_COLORS = {"CCEMG": "#D55E00", ABC_METHOD: "#0072B2"}


@dataclass
class CCEFit:
    betas: np.ndarray                 # N x k
    x_tilde: np.ndarray               # T x N x k
    residuals: np.ndarray             # T x N
    gram: np.ndarray                  # N x k x k, Xi'MXi/T
    intercept_y: np.ndarray           # N
    proxy_y: np.ndarray               # N x (k+1)
    intercept_x: np.ndarray           # N x k
    proxy_x: np.ndarray               # N x (k+1) x k
    residual_se: np.ndarray           # N
    y_bar: np.ndarray                 # T
    x_bar: np.ndarray                 # T x k
    rank_joint: np.ndarray            # N


@dataclass
class ABCFit:
    beta_abc: np.ndarray
    b0_hat: np.ndarray
    b1_hat: np.ndarray
    b0_unit: np.ndarray
    b1_unit: np.ndarray
    unit_corrected_slopes: np.ndarray
    pbar: np.ndarray


def fit_cce(y: np.ndarray, x: np.ndarray, rows: np.ndarray | None = None) -> CCEFit:
    """Fit individual CCE regressions for one state.

    y is T x N and x is T x N x k.  The augmentation contains an intercept,
    the cross-sectional average of y, and all k cross-sectional regressor
    averages, as in Section 5.
    """
    if rows is not None:
        y = y[rows]
        x = x[rows]
    t_count, n_count = y.shape
    k_count = x.shape[2]
    y_bar = y.mean(axis=1)
    x_bar = x.mean(axis=1)
    controls = np.column_stack((np.ones(t_count), y_bar, x_bar))

    betas = np.empty((n_count, k_count))
    x_tilde = np.empty((t_count, n_count, k_count))
    residuals = np.empty((t_count, n_count))
    gram = np.empty((n_count, k_count, k_count))
    intercept_y = np.empty(n_count)
    proxy_y = np.empty((n_count, k_count + 1))
    intercept_x = np.empty((n_count, k_count))
    proxy_x = np.empty((n_count, k_count + 1, k_count))
    residual_se = np.empty(n_count)
    rank_joint = np.empty(n_count, dtype=int)

    for index in range(n_count):
        x_i = x[:, index, :]
        y_i = y[:, index]
        coef_x = np.linalg.lstsq(controls, x_i, rcond=1e-10)[0]
        coef_y_controls = np.linalg.lstsq(controls, y_i, rcond=1e-10)[0]
        x_res = x_i - controls @ coef_x
        y_res = y_i - controls @ coef_y_controls
        beta_i = np.linalg.pinv(x_res.T @ x_res, rcond=1e-10) @ (x_res.T @ y_res)
        y_nuisance = np.linalg.lstsq(controls, y_i - x_i @ beta_i, rcond=1e-10)[0]
        error_i = y_i - x_i @ beta_i - controls @ y_nuisance
        joint = np.column_stack((controls, x_i))
        rank_i = int(np.linalg.matrix_rank(joint, tol=1e-10))
        dof = max(1, t_count - rank_i)

        betas[index] = beta_i
        x_tilde[:, index, :] = x_res
        residuals[:, index] = error_i
        gram[index] = x_res.T @ x_res / t_count
        intercept_y[index] = y_nuisance[0]
        proxy_y[index] = y_nuisance[1:]
        intercept_x[index] = coef_x[0]
        proxy_x[index] = coef_x[1:]
        residual_se[index] = np.sqrt(np.sum(error_i**2) / dof)
        rank_joint[index] = rank_i

    return CCEFit(
        betas=betas,
        x_tilde=x_tilde,
        residuals=residuals,
        gram=gram,
        intercept_y=intercept_y,
        proxy_y=proxy_y,
        intercept_x=intercept_x,
        proxy_x=proxy_x,
        residual_se=residual_se,
        y_bar=y_bar,
        x_bar=x_bar,
        rank_joint=rank_joint,
    )


def analytical_bias_correction(fit: CCEFit, y: np.ndarray, x: np.ndarray,
                               factor_count: int = 2) -> ABCFit:
    """Compute the vector ABC estimator with the subtraction sign in Eq. (44).

    The simulation design has two common factors.  With seven empirical
    regressors, the full CCE augmentation has eight averages, so the nuisance
    loading space used in the analytical correction is reduced to the requested
    factor_count using principal components of (y_bar, x_bar).  The CCEMG
    slopes themselves continue to use all cross-sectional averages.
    """
    n_count, k_count = fit.betas.shape
    t_count = fit.residuals.shape[0]
    augmented_dimension = k_count + 2
    if factor_count < 1 or factor_count > min(k_count + 1, t_count - 1):
        raise ValueError("factor_count is outside the admissible range.")

    averages = np.column_stack((fit.y_bar, fit.x_bar))
    averages = averages - averages.mean(axis=0)
    left_vectors = np.linalg.svd(averages, full_matrices=False)[0]
    factors = left_vectors[:, :factor_count] * np.sqrt(t_count)
    factor_controls = np.column_stack((np.ones(t_count), factors))

    p_units = np.empty((n_count, factor_count + 1, augmented_dimension))
    pi_units = np.empty((n_count, factor_count + 1))
    Pi_units = np.empty((n_count, factor_count + 1, k_count))

    for index in range(n_count):
        beta_i = fit.betas[index]
        coefficients_x = np.linalg.lstsq(factor_controls, x[:, index, :], rcond=1e-10)[0]
        coefficients_y = np.linalg.lstsq(
            factor_controls, y[:, index] - x[:, index, :] @ beta_i, rcond=1e-10
        )[0]
        A_i = coefficients_x[0]
        Gamma_i = coefficients_x[1:]
        a_i = coefficients_y[0]
        gamma_i = coefficients_y[1:]
        B_i = np.concatenate(([a_i + A_i @ beta_i], A_i))
        C_i = np.column_stack((gamma_i + Gamma_i @ beta_i, Gamma_i))
        p_units[index] = np.block([
            [np.concatenate(([1.0], B_i))[None, :]],
            [np.column_stack((np.zeros(factor_count), C_i))],
        ])
        pi_units[index] = np.concatenate(([a_i], gamma_i))
        Pi_units[index] = np.vstack((A_i, Gamma_i))

    pbar = p_units.mean(axis=0)
    pbar_plus = np.linalg.pinv(pbar, rcond=1e-10)
    b0_unit = np.empty((n_count, k_count))
    b1_unit = np.empty((n_count, k_count))

    for index in range(n_count):
        beta_i = fit.betas[index]
        selector = np.column_stack((
            np.zeros(k_count),
            beta_i,
            np.eye(k_count),
        ))
        sigma_vector = np.concatenate((
            [0.0, fit.residual_se[index]],
            np.ones(k_count),
        ))
        inverse_gram = np.linalg.pinv(fit.gram[index], rcond=1e-10)
        b0_unit[index] = selector @ pbar_plus @ pi_units[index]
        b1_unit[index] = (
            inverse_gram
            @ Pi_units[index].T
            @ pbar_plus.T
            @ sigma_vector
        )

    b0_hat = np.sqrt(t_count / n_count) * b0_unit.mean(axis=0)
    b1_hat = np.sqrt(t_count / n_count) * b1_unit.mean(axis=0)
    beta_mg = fit.betas.mean(axis=0)
    beta_abc = beta_mg - (b0_hat + b1_hat) / np.sqrt(n_count * t_count)

    unit_corrected_slopes = fit.betas - (b0_unit + b1_unit) / n_count
    equivalent = unit_corrected_slopes.mean(axis=0)
    if not np.allclose(beta_abc, equivalent, rtol=1e-9, atol=1e-9):
        raise AssertionError("ABC scaling identity failed.")
    return ABCFit(
        beta_abc, b0_hat, b1_hat, b0_unit, b1_unit,
        unit_corrected_slopes, pbar,
    )


def standard_ccemg_covariance(betas: np.ndarray) -> np.ndarray:
    n_count = betas.shape[0]
    centered = betas - betas.mean(axis=0)
    return centered.T @ centered / (n_count * (n_count - 1))


def whole_sum_dwb_l1(
    fit: CCEFit,
    abc: ABCFit,
    bootstrap_reps: int,
    confidence: float,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Whole-sum cross-sectional multiplier bootstrap with l_N fixed at one.

    For vector slopes, row i of unit_influence is
    sqrt(N/(N-1))*sqrt(T)*(beta_ABC,i-beta_ABC). The conditional covariance of
    N^(-1/2) sum_i unit_influence_i xi_i is therefore Zc'Zc/N. Dividing by NT
    gives the covariance of beta_ABC.
    """
    t_count, n_count = fit.residuals.shape
    if n_count <= 1:
        raise ValueError("The l_N=1 covariance requires at least two units.")

    centered = abc.unit_corrected_slopes - abc.beta_abc
    unit_influence = np.sqrt(t_count) * np.sqrt(n_count / (n_count - 1)) * centered
    multipliers = rng.normal(size=(n_count, bootstrap_reps))
    statistic_draws = multipliers.T @ unit_influence / np.sqrt(n_count)

    conditional_covariance = unit_influence.T @ unit_influence / n_count
    estimator_covariance = conditional_covariance / (n_count * t_count)
    corrected_slope_covariance = centered.T @ centered / (n_count * (n_count - 1))
    if not np.allclose(estimator_covariance, corrected_slope_covariance,
                       rtol=1e-10, atol=1e-12):
        raise AssertionError("The l_N=1 covariance identity failed.")
    bootstrap_covariance = np.cov(statistic_draws, rowvar=False, ddof=1) / (
        n_count * t_count
    )

    alpha = 1 - confidence
    q_low = np.quantile(statistic_draws, alpha / 2, axis=0, method="median_unbiased")
    q_high = np.quantile(statistic_draws, 1 - alpha / 2, axis=0, method="median_unbiased")
    scale = np.sqrt(n_count * t_count)
    lower = abc.beta_abc - q_high / scale
    upper = abc.beta_abc - q_low / scale
    se = np.sqrt(np.maximum(np.diag(estimator_covariance), 0))
    return {
        "lower": lower,
        "upper": upper,
        "se": se,
        "covariance": estimator_covariance,
        "conditional_covariance": conditional_covariance,
        "bootstrap_covariance": bootstrap_covariance,
        "q_low": q_low,
        "q_high": q_high,
        "unit_influence": unit_influence,
        "bandwidth": 1,
        "finite_sample_correction": True,
        "covariance_identity_max_abs_error": float(
            np.max(np.abs(estimator_covariance - corrected_slope_covariance))
        ),
        "bootstrap_covariance_max_abs_difference": float(
            np.max(np.abs(bootstrap_covariance - estimator_covariance))
        ),
    }


def load_state_panel(data: pd.DataFrame, state: str) -> tuple[np.ndarray, np.ndarray, list[int], list[str]]:
    subset = data.loc[data["State_norm"] == state].copy()
    years = sorted(subset["Year"].unique().tolist())
    units = sorted(subset["Id"].unique().tolist())
    expected = pd.MultiIndex.from_product((units, years), names=("Id", "Year"))
    indexed = subset.set_index(["Id", "Year"]).sort_index()
    if indexed.index.has_duplicates:
        raise ValueError(f"{state}: duplicate Id-Year observations.")
    if not indexed.index.equals(expected):
        missing = expected.difference(indexed.index)
        raise ValueError(f"{state}: panel is not balanced; missing {len(missing)} cells.")
    if indexed[["ln_yields", *VARIABLES]].isna().any().any():
        raise ValueError(f"{state}: missing values in required model variables.")

    t_count, n_count, k_count = len(years), len(units), len(VARIABLES)
    y = np.empty((t_count, n_count))
    x = np.empty((t_count, n_count, k_count))
    region_names: list[str] = []
    for unit_index, unit in enumerate(units):
        unit_data = indexed.loc[unit].reindex(years)
        y[:, unit_index] = unit_data["ln_yields"].to_numpy(float)
        x[:, unit_index, :] = unit_data[VARIABLES].to_numpy(float)
        region_names.append(str(unit_data["Region"].iloc[0]))
    return y, x, years, region_names


def estimate_state(
    data: pd.DataFrame,
    state: str,
    bootstrap_reps: int,
    confidence: float,
    factor_count: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    y, x, years, region_names = load_state_panel(data, state)
    t_count, n_count = y.shape
    fit = fit_cce(y, x)
    beta_mg = fit.betas.mean(axis=0)
    abc = analytical_bias_correction(fit, y, x, factor_count)
    standard_cov = standard_ccemg_covariance(fit.betas)
    standard_se = np.sqrt(np.maximum(np.diag(standard_cov), 0))
    z_value = norm.ppf((1 + confidence) / 2)
    standard_lower = beta_mg - z_value * standard_se
    standard_upper = beta_mg + z_value * standard_se

    bootstrap = whole_sum_dwb_l1(
        fit, abc, bootstrap_reps, confidence, np.random.default_rng(seed)
    )

    rows: list[dict[str, object]] = []
    for variable_index, variable in enumerate(VARIABLES):
        for method, estimate, se, lower, upper in (
            ("CCEMG", beta_mg[variable_index], standard_se[variable_index],
             standard_lower[variable_index], standard_upper[variable_index]),
            (ABC_METHOD, abc.beta_abc[variable_index], bootstrap["se"][variable_index],
             bootstrap["lower"][variable_index], bootstrap["upper"][variable_index]),
        ):
            rows.append({
                "State": state,
                "Variable": variable,
                "Variable_label": DISPLAY_NAMES[variable],
                "Method": method,
                "Estimate": estimate,
                "Standard_error": se,
                "CI_lower": lower,
                "CI_upper": upper,
                "CI_length": upper - lower,
                "p_value_normal": 2 * norm.sf(abs(estimate / se)) if se > 0 else np.nan,
                "N": n_count,
                "T": t_count,
                "Confidence_level": confidence,
                "Bootstrap_reps": bootstrap_reps if method != "CCEMG" else np.nan,
                "DWB_bandwidth": bootstrap["bandwidth"] if method != "CCEMG" else np.nan,
            })
    results = pd.DataFrame(rows)

    diagnostics = pd.DataFrame([{
        "State": state,
        "N": n_count,
        "T": t_count,
        "Start_year": min(years),
        "End_year": max(years),
        "Regions": " | ".join(region_names),
        "Joint_rank_min": int(fit.rank_joint.min()),
        "Joint_rank_max": int(fit.rank_joint.max()),
        "Pbar_rank": int(np.linalg.matrix_rank(abc.pbar, tol=1e-10)),
        "Pbar_condition_number": float(np.linalg.cond(abc.pbar)),
        "ABC_factor_count": int(factor_count),
        "DWB_bandwidth": int(bootstrap["bandwidth"]),
        "DWB_dimension": "cross_section",
        "DWB_finite_sample_correction": bool(bootstrap["finite_sample_correction"]),
        "Conditional_covariance_min_eigenvalue": float(
            np.linalg.eigvalsh(bootstrap["conditional_covariance"]).min()
        ),
        "Conditional_covariance_max_eigenvalue": float(
            np.linalg.eigvalsh(bootstrap["conditional_covariance"]).max()
        ),
        "Covariance_identity_max_abs_error": bootstrap["covariance_identity_max_abs_error"],
        "Bootstrap_covariance_max_abs_difference": bootstrap[
            "bootstrap_covariance_max_abs_difference"
        ],
    }])

    slope_rows = []
    for unit_index, region in enumerate(region_names):
        for variable_index, variable in enumerate(VARIABLES):
            slope_rows.append({
                "State": state,
                "Region": region,
                "Variable": variable,
                "Individual_CCE_slope": fit.betas[unit_index, variable_index],
                "ABC_b0_unit_component": abc.b0_unit[unit_index, variable_index],
                "ABC_b1_unit_component": abc.b1_unit[unit_index, variable_index],
                "Individual_ABC_corrected_slope": abc.unit_corrected_slopes[
                    unit_index, variable_index
                ],
                "DWB_centered_contribution": bootstrap["unit_influence"][
                    unit_index, variable_index
                ],
            })
    return results, diagnostics, pd.DataFrame(slope_rows)


def plot_state(results: pd.DataFrame, state: str, destination: Path) -> None:
    state_data = results.loc[results["State"] == state].copy()
    figure, axes = plt.subplots(len(VARIABLES), 1, figsize=(10.5, 12.5), constrained_layout=True)
    methods = ["CCEMG", ABC_METHOD]
    y_positions = {"CCEMG": 1.0, ABC_METHOD: 0.0}
    for axis, variable in zip(axes, VARIABLES):
        variable_data = state_data.loc[state_data["Variable"] == variable]
        for method in methods:
            row = variable_data.loc[variable_data["Method"] == method].iloc[0]
            y_position = y_positions[method]
            color = METHOD_COLORS[method]
            axis.fill_betweenx(
                [y_position - 0.18, y_position + 0.18],
                row["CI_lower"], row["CI_upper"],
                color=color, alpha=0.22, linewidth=0,
            )
            axis.plot(row["Estimate"], y_position, marker="o", color=color, markersize=7)
        axis.axvline(0, color="#666666", linestyle="--", linewidth=0.9)
        axis.set_yticks([0, 1], ["ABC + DWB, l_N=1", "CCEMG"])
        axis.set_title(DISPLAY_NAMES[variable], loc="left", fontsize=10.5, fontweight="semibold")
        axis.grid(axis="x", color="#dddddd", linewidth=0.7)
        axis.set_axisbelow(True)
        for spine in ("top", "right", "left"):
            axis.spines[spine].set_visible(False)
    axes[-1].set_xlabel("Estimated mean effect with 95% confidence interval")
    figure.suptitle(f"{state}: CCEMG and ABC whole-sum DWB mean effects (l_N=1)", fontsize=15)
    figure.savefig(destination, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def write_summary(results: pd.DataFrame, diagnostics: pd.DataFrame, destination: Path) -> None:
    lines = [
        "# Empirical CCEMG and ABC whole-sum DWB results (l_N = 1)",
        "",
        "The estimates use the Section 5 wheat-yield specification separately for NSW, SA, and VIC.",
        "CCEMG intervals use Pesaran's cross-sectional dispersion estimator. ABC standard errors",
        "use the exact finite-sample-corrected l_N=1 covariance of the corrected individual slopes.",
        "ABC confidence intervals use asymmetric whole-sum multiplier-bootstrap quantiles, with",
        "the analytical correction subtracted as in Equation (44). JCK is not estimated.",
        "",
    ]
    for state in STATES:
        lines.extend([f"## {state}", "", "| Variable | Method | Estimate | SE | CI lower | CI upper | CI length |", "|---|---|---:|---:|---:|---:|---:|"])
        for _, row in results.loc[results["State"] == state].iterrows():
            lines.append(
                f"| {row['Variable_label']} | {row['Method']} | {row['Estimate']:.6f} | "
                f"{row['Standard_error']:.6f} | {row['CI_lower']:.6f} | "
                f"{row['CI_upper']:.6f} | {row['CI_length']:.6f} |"
            )
        lines.append("")
    diagnostic_columns = [
        "State", "N", "T", "Start_year", "End_year", "Joint_rank_min",
        "Joint_rank_max", "Pbar_rank", "Pbar_condition_number", "DWB_bandwidth",
        "DWB_dimension", "DWB_finite_sample_correction",
        "Covariance_identity_max_abs_error", "Bootstrap_covariance_max_abs_difference",
        "ABC_factor_count",
    ]
    lines.extend([
        "## Diagnostics", "",
        "| " + " | ".join(diagnostic_columns) + " |",
        "|" + "|".join(["---"] * len(diagnostic_columns)) + "|",
    ])
    for _, row in diagnostics.iterrows():
        formatted = []
        for column in diagnostic_columns:
            value = row[column]
            formatted.append(f"{value:.6g}" if isinstance(value, (float, np.floating)) else str(value))
        lines.append("| " + " | ".join(formatted) + " |")
    lines.append("")
    destination.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    default_data = Path(__file__).resolve().parent / "data" / "Aus_data_with_agg_new.xlsx"
    parser.add_argument("--data", type=Path, default=default_data)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--bootstrap-reps", type=int, default=9999)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--factor-count", type=int, default=2,
                        help="Number of common factors used in the ABC nuisance correction.")
    parser.add_argument("--seed", type=int, default=20261001)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.bootstrap_reps < 99:
        raise ValueError("Use at least 99 bootstrap replications.")
    if not 0 < args.confidence < 1:
        raise ValueError("confidence must lie in (0,1).")
    data = pd.read_excel(args.data)
    required = {"Id", "Region", "Year", "State", "ln_yields", *VARIABLES}
    missing = required.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    data["State_norm"] = data["State"].astype(str).str.strip().str.upper()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_results, all_diagnostics, all_slopes = [], [], []
    for state_index, state in enumerate(STATES):
        results, diagnostics, slopes = estimate_state(
            data=data,
            state=state,
            bootstrap_reps=args.bootstrap_reps,
            confidence=args.confidence,
            factor_count=args.factor_count,
            seed=args.seed + 1000 * state_index,
        )
        all_results.append(results)
        all_diagnostics.append(diagnostics)
        all_slopes.append(slopes)

    results = pd.concat(all_results, ignore_index=True)
    diagnostics = pd.concat(all_diagnostics, ignore_index=True)
    slopes = pd.concat(all_slopes, ignore_index=True)
    results.to_csv(args.output_dir / "empirical_ccemg_abc_results.csv", index=False)
    diagnostics.to_csv(args.output_dir / "empirical_diagnostics.csv", index=False)
    slopes.to_csv(args.output_dir / "individual_cce_slopes_and_abc_components.csv", index=False)
    lengths = results.pivot_table(
        index=["State", "Variable", "Variable_label"],
        columns="Method", values="CI_length",
    ).reset_index()
    lengths.to_csv(args.output_dir / "confidence_interval_lengths.csv", index=False)
    for state in STATES:
        plot_state(results, state, args.output_dir / f"{state}_mean_effects_with_shaded_CI.png")
    write_summary(results, diagnostics, args.output_dir / "RESULTS_SUMMARY.md")
    run_info = {
        "data_file": str(args.data),
        "states": list(STATES),
        "variables": VARIABLES,
        "bootstrap_reps": args.bootstrap_reps,
        "confidence": args.confidence,
        "dwb_bandwidth": 1,
        "dwb_dimension": "cross_section",
        "finite_sample_correction": True,
        "factor_count": args.factor_count,
        "seed": args.seed,
        "interval": "asymmetric basic whole-sum multiplier bootstrap",
        "abc_variance": "exact l_N=1 finite-sample-corrected covariance",
        "abc_sign": "subtraction (Equation 44)",
    }
    (args.output_dir / "run_configuration.json").write_text(
        json.dumps(run_info, indent=2), encoding="utf-8"
    )
    print(results.to_string(index=False))
    print(f"\nWrote results to {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
