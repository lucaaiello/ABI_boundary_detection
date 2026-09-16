from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import ks_2samp, wasserstein_distance


PARAMETERS = ["beta0", "sigma2_w", "eta_raw", "eta", "rho"]
JOINT_PARAMETERS = ["beta0", "sigma2_w", "eta", "rho"]
PARAMETER_LABELS = {
    "beta0": r"$\beta_0$",
    "sigma2_w": r"$\sigma_w^2$",
    "eta_raw": r"$\eta_{\mathrm{raw}}$",
    "eta": r"$\eta$",
    "rho": r"$\rho$",
}
ABI_COLOR = "#2f6f8f"
MCMC_COLOR = "#d77a3d"
REFERENCE_COLOR = "#252525"
INTERACTION_COLOR = "#67805f"


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")
    return path


def interval_overlap_fraction(row: pd.Series) -> float:
    overlap = max(
        0.0,
        min(row["upper_95_abi"], row["upper_95_mcmc"])
        - max(row["lower_95_abi"], row["lower_95_mcmc"]),
    )
    smaller_width = min(
        row["upper_95_abi"] - row["lower_95_abi"],
        row["upper_95_mcmc"] - row["lower_95_mcmc"],
    )
    return overlap / smaller_width if smaller_width > 0 else np.nan


def load_parameter_summaries(results_dir: Path, method: str) -> pd.DataFrame:
    frame = pd.read_csv(require_file(results_dir / "combined_parameter_summaries.csv"))
    frame["method"] = method
    return frame


def summarize_methods(parameters: pd.DataFrame) -> pd.DataFrame:
    frame = parameters.copy()
    frame["interval_width"] = frame["upper_95"] - frame["lower_95"]
    return (
        frame.groupby(["method", "parameter"], as_index=False)
        .agg(
            n_datasets=("dataset_id", "nunique"),
            mean_bias=("bias_mean", "mean"),
            mean_abs_error=("abs_error_mean", "mean"),
            median_abs_error=("abs_error_mean", "median"),
            coverage_95=("covered_95", "mean"),
            mean_posterior_sd=("posterior_sd", "mean"),
            mean_interval_width=("interval_width", "mean"),
        )
    )


def load_boundary_metric_summaries(
    abi_results_dir: Path,
    mcmc_results_dir: Path,
    dataset_ids: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    for results_dir, method in (
        (abi_results_dir, "ABI"),
        (mcmc_results_dir, "Matched MCMC"),
    ):
        frame = pd.read_csv(require_file(results_dir / "combined_edge_metrics.csv"))
        frame = frame.loc[frame["dataset_id"].isin(dataset_ids)].copy()
        frame["method"] = method
        frames.append(frame)

    by_dataset = pd.concat(frames, ignore_index=True)
    metrics = [
        "auroc",
        "average_precision",
        "brier",
        "sensitivity_mpm",
        "specificity_mpm",
        "posterior_boundary_count_mpm",
        "true_boundary_count",
        "boundary_count_mean_draws",
        "boundary_count_truth_in_95",
    ]
    summary_rows: list[dict[str, float | int | str]] = []
    for method, group in by_dataset.groupby("method", sort=False):
        row: dict[str, float | int | str] = {
            "method": method,
            "n_datasets": int(group["dataset_id"].nunique()),
        }
        for metric in metrics:
            values = group[metric].dropna()
            row[f"{metric}_n"] = int(len(values))
            row[f"{metric}_mean"] = float(values.mean()) if len(values) else np.nan
            row[f"{metric}_median"] = float(values.median()) if len(values) else np.nan
        summary_rows.append(row)
    return by_dataset, pd.DataFrame(summary_rows)


def pair_parameter_summaries(abi: pd.DataFrame, mcmc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    columns = [
        "dataset_id",
        "parameter",
        "posterior_mean",
        "posterior_sd",
        "posterior_median",
        "lower_95",
        "upper_95",
        "truth",
    ]
    paired = abi[columns].merge(
        mcmc[columns], on=["dataset_id", "parameter"], suffixes=("_abi", "_mcmc"), validate="one_to_one"
    )
    if not np.allclose(paired["truth_abi"], paired["truth_mcmc"], rtol=0.0, atol=1e-10):
        raise ValueError("ABI and matched MCMC summaries contain different parameter truths.")
    paired["posterior_mean_difference"] = paired["posterior_mean_abi"] - paired["posterior_mean_mcmc"]
    paired["absolute_posterior_mean_difference"] = paired["posterior_mean_difference"].abs()
    paired["posterior_median_difference"] = paired["posterior_median_abi"] - paired["posterior_median_mcmc"]
    paired["absolute_posterior_median_difference"] = paired["posterior_median_difference"].abs()
    paired["posterior_width_ratio_abi_to_mcmc"] = (
        (paired["upper_95_abi"] - paired["lower_95_abi"])
        / (paired["upper_95_mcmc"] - paired["lower_95_mcmc"]).replace(0.0, np.nan)
    )
    paired["interval_overlap_fraction"] = paired.apply(interval_overlap_fraction, axis=1)
    summary = (
        paired.groupby("parameter", as_index=False)
        .agg(
            n_datasets=("dataset_id", "nunique"),
            mean_abs_posterior_mean_difference=("absolute_posterior_mean_difference", "mean"),
            median_abs_posterior_mean_difference=("absolute_posterior_mean_difference", "median"),
            mean_abs_posterior_median_difference=("absolute_posterior_median_difference", "mean"),
            median_width_ratio_abi_to_mcmc=("posterior_width_ratio_abi_to_mcmc", "median"),
            mean_interval_overlap_fraction=("interval_overlap_fraction", "mean"),
        )
    )
    correlations = (
        paired.groupby("parameter")
        .apply(
            lambda group: group["posterior_mean_abi"].corr(group["posterior_mean_mcmc"]),
            include_groups=False,
        )
        .rename("posterior_mean_correlation")
        .reset_index()
    )
    return paired, summary.merge(correlations, on="parameter", how="left")


def posterior_draw_path(results_dir: Path, dataset_id: str) -> Path:
    return (
        results_dir
        / "per_dataset"
        / dataset_id
        / f"{dataset_id}_posterior_draws.csv.gz"
    )


def compare_posterior_draws(
    abi_results_dir: Path, mcmc_results_dir: Path, dataset_ids: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    marginal_rows: list[dict[str, float | str]] = []
    joint_rows: list[dict[str, float | str]] = []
    for index, dataset_id in enumerate(dataset_ids, start=1):
        abi_draws = pd.read_csv(require_file(posterior_draw_path(abi_results_dir, dataset_id)))
        mcmc_draws = pd.read_csv(require_file(posterior_draw_path(mcmc_results_dir, dataset_id)))
        for parameter in PARAMETERS:
            abi_values = abi_draws[parameter].to_numpy(dtype=float)
            mcmc_values = mcmc_draws[parameter].to_numpy(dtype=float)
            mcmc_sd = float(np.std(mcmc_values, ddof=1))
            distance = float(wasserstein_distance(abi_values, mcmc_values))
            marginal_rows.append(
                {
                    "dataset_id": dataset_id,
                    "parameter": parameter,
                    "wasserstein_distance": distance,
                    "wasserstein_over_mcmc_sd": distance / mcmc_sd if mcmc_sd > 0 else np.nan,
                    "ks_statistic": float(ks_2samp(abi_values, mcmc_values).statistic),
                }
            )

        abi_correlation = abi_draws[JOINT_PARAMETERS].corr().to_numpy()
        mcmc_correlation = mcmc_draws[JOINT_PARAMETERS].corr().to_numpy()
        difference = abi_correlation - mcmc_correlation
        joint_rows.append(
            {
                "dataset_id": dataset_id,
                "correlation_frobenius_difference": float(np.linalg.norm(difference, ord="fro")),
                "maximum_absolute_correlation_difference": float(np.max(np.abs(difference))),
            }
        )
        if index % 10 == 0 or index == len(dataset_ids):
            print(f"Compared posterior draws for {index}/{len(dataset_ids)} datasets")

    marginal = pd.DataFrame(marginal_rows)
    marginal_summary = (
        marginal.groupby("parameter", as_index=False)
        .agg(
            n_datasets=("dataset_id", "nunique"),
            mean_wasserstein=("wasserstein_distance", "mean"),
            median_wasserstein=("wasserstein_distance", "median"),
            mean_scaled_wasserstein=("wasserstein_over_mcmc_sd", "mean"),
            median_scaled_wasserstein=("wasserstein_over_mcmc_sd", "median"),
            mean_ks=("ks_statistic", "mean"),
            median_ks=("ks_statistic", "median"),
        )
    )
    return marginal, marginal_summary, pd.DataFrame(joint_rows)


def load_edge_pairs(
    abi_results_dir: Path, mcmc_results_dir: Path, dataset_ids: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pair_frames: list[pd.DataFrame] = []
    summary_rows: list[dict[str, float | int | str]] = []
    for dataset_id in dataset_ids:
        abi_path = (
            abi_results_dir / "per_dataset" / dataset_id / f"{dataset_id}_edge_probabilities.csv"
        )
        mcmc_path = (
            mcmc_results_dir / "per_dataset" / dataset_id / f"{dataset_id}_edge_probabilities.csv"
        )
        abi = pd.read_csv(require_file(abi_path))
        mcmc = pd.read_csv(require_file(mcmc_path))
        edge_keys = [
            "edge_index_1based",
            "node_i_1based",
            "node_j_1based",
            "boundary_true",
        ]
        paired = abi.merge(
            mcmc[edge_keys + ["boundary_prob_mcmc", "boundary_median_mcmc"]],
            on=edge_keys,
            validate="one_to_one",
        )
        if len(paired) != len(abi) or len(paired) != len(mcmc):
            raise ValueError(f"ABI and matched MCMC edge identities differ for {dataset_id}.")
        probability_difference = paired["boundary_prob_abi"] - paired["boundary_prob_mcmc"]
        abi_selected = paired["boundary_median_abi"].astype(bool)
        mcmc_selected = paired["boundary_median_mcmc"].astype(bool)
        union = int((abi_selected | mcmc_selected).sum())
        shared = int((abi_selected & mcmc_selected).sum())
        summary_rows.append(
            {
                "dataset_id": dataset_id,
                "edge_count": len(paired),
                "boundary_probability_correlation": paired["boundary_prob_abi"].corr(
                    paired["boundary_prob_mcmc"]
                ),
                "boundary_probability_mae": float(probability_difference.abs().mean()),
                "boundary_probability_rmse": float(np.sqrt(np.mean(probability_difference**2))),
                "abi_selected": int(abi_selected.sum()),
                "mcmc_selected": int(mcmc_selected.sum()),
                "shared_selected": shared,
                "jaccard_selected": shared / union if union > 0 else 1.0,
            }
        )
        pair_frames.append(paired)
    return pd.concat(pair_frames, ignore_index=True), pd.DataFrame(summary_rows)


def build_error_decomposition(
    parameter_pairs: pd.DataFrame,
    diagnostics: pd.DataFrame,
    abi_runtime: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    precision = diagnostics[
        ["dataset_id", "parameter", "effective_sample_size", "mcse_mean"]
    ]
    abi_draw_counts = abi_runtime[["dataset_id", "n_saved_draws"]].rename(
        columns={"n_saved_draws": "n_abi_draws"}
    )
    frame = parameter_pairs.loc[
        parameter_pairs["parameter"].isin(JOINT_PARAMETERS)
    ].merge(precision, on=["dataset_id", "parameter"], validate="one_to_one").merge(
        abi_draw_counts, on="dataset_id", validate="many_to_one"
    )
    frame["abi_error"] = frame["posterior_mean_abi"] - frame["truth_abi"]
    frame["mcmc_error"] = frame["posterior_mean_mcmc"] - frame["truth_mcmc"]
    frame["approximation_error"] = (
        frame["posterior_mean_abi"] - frame["posterior_mean_mcmc"]
    )
    frame["abi_squared_error"] = frame["abi_error"] ** 2
    frame["mcmc_squared_error"] = frame["mcmc_error"] ** 2
    frame["approximation_squared_error"] = frame["approximation_error"] ** 2
    frame["interaction"] = 2.0 * frame["mcmc_error"] * frame["approximation_error"]
    frame["mcmc_mc_variance"] = frame["mcse_mean"] ** 2
    frame["abi_mc_variance"] = frame["posterior_sd_abi"] ** 2 / frame["n_abi_draws"]

    summary_rows: list[dict[str, float | int | str]] = []
    for parameter in JOINT_PARAMETERS:
        group = frame.loc[frame["parameter"] == parameter]
        if group.empty:
            continue
        abi_mse = float(group["abi_squared_error"].mean())
        mcmc_mse = float(group["mcmc_squared_error"].mean())
        approximation_mse = float(group["approximation_squared_error"].mean())
        interaction = float(group["interaction"].mean())
        mcmc_mc_variance = float(group["mcmc_mc_variance"].mean())
        abi_mc_variance = float(group["abi_mc_variance"].mean())
        corrected_abi_mse = abi_mse - abi_mc_variance
        corrected_mcmc_mse = mcmc_mse - mcmc_mc_variance
        corrected_approximation_mse = (
            approximation_mse - mcmc_mc_variance - abi_mc_variance
        )
        corrected_interaction = interaction + 2.0 * mcmc_mc_variance
        reconstructed = (
            corrected_mcmc_mse
            + corrected_approximation_mse
            + corrected_interaction
        )
        summary_rows.append(
            {
                "parameter": parameter,
                "n_datasets": int(group["dataset_id"].nunique()),
                "abi_total_mse_observed": abi_mse,
                "abi_total_mse_mc_corrected": corrected_abi_mse,
                "mcmc_reference_mse_observed": mcmc_mse,
                "approximation_mse_observed": approximation_mse,
                "interaction_observed": interaction,
                "mean_abi_mc_variance": abi_mc_variance,
                "mean_mcmc_mc_variance": mcmc_mc_variance,
                "mcmc_reference_mse_mc_corrected": corrected_mcmc_mse,
                "approximation_mse_mc_corrected": corrected_approximation_mse,
                "interaction_mc_corrected": corrected_interaction,
                "reconstructed_abi_mse": reconstructed,
                "decomposition_residual": corrected_abi_mse - reconstructed,
                "mc_variance_pct_of_observed_approximation_mse": (
                    100.0 * (abi_mc_variance + mcmc_mc_variance) / approximation_mse
                    if approximation_mse > 0.0
                    else np.nan
                ),
                "mcse_of_aggregate_abi_bias": float(
                    np.sqrt(group["abi_mc_variance"].sum()) / len(group)
                ),
                "mcse_of_aggregate_mcmc_bias": float(
                    np.sqrt(group["mcmc_mc_variance"].sum()) / len(group)
                ),
            }
        )
    return frame, pd.DataFrame(summary_rows)


def run_comparison(
    benchmark_dir: Path,
    abi_results_dir: Path,
    mcmc_results_dir: Path,
    output_dir: Path,
) -> dict[str, pd.DataFrame]:
    benchmark_dir = Path(benchmark_dir).resolve()
    abi_results_dir = Path(abi_results_dir).resolve()
    mcmc_results_dir = Path(mcmc_results_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(require_file(benchmark_dir / "benchmark_manifest.csv"))
    abi_parameters = load_parameter_summaries(abi_results_dir, "ABI")
    mcmc_parameters = load_parameter_summaries(mcmc_results_dir, "Matched MCMC")
    common_ids = sorted(
        set(manifest["dataset_id"])
        & set(abi_parameters["dataset_id"])
        & set(mcmc_parameters["dataset_id"])
    )
    if not common_ids:
        raise ValueError("ABI and matched MCMC have no datasets in common.")

    abi_parameters = abi_parameters[abi_parameters["dataset_id"].isin(common_ids)].copy()
    mcmc_parameters = mcmc_parameters[mcmc_parameters["dataset_id"].isin(common_ids)].copy()
    parameter_all = pd.concat([abi_parameters, mcmc_parameters], ignore_index=True)
    method_summary = summarize_methods(parameter_all)
    parameter_pairs, parameter_agreement = pair_parameter_summaries(
        abi_parameters, mcmc_parameters
    )
    marginal_distances, marginal_summary, joint_differences = compare_posterior_draws(
        abi_results_dir, mcmc_results_dir, common_ids
    )
    edge_pairs, edge_agreement = load_edge_pairs(abi_results_dir, mcmc_results_dir, common_ids)
    edge_metrics, edge_metric_summary = load_boundary_metric_summaries(
        abi_results_dir, mcmc_results_dir, common_ids
    )

    edge_agreement_summary = pd.DataFrame(
        [
            {
                "n_datasets": len(edge_agreement),
                "mean_probability_correlation": edge_agreement[
                    "boundary_probability_correlation"
                ].mean(),
                "median_probability_correlation": edge_agreement[
                    "boundary_probability_correlation"
                ].median(),
                "mean_probability_mae": edge_agreement["boundary_probability_mae"].mean(),
                "median_probability_mae": edge_agreement["boundary_probability_mae"].median(),
                "mean_jaccard": edge_agreement["jaccard_selected"].mean(),
                "median_jaccard": edge_agreement["jaccard_selected"].median(),
            }
        ]
    )
    abi_runtime = pd.read_csv(require_file(abi_results_dir / "combined_runtime.csv"))
    mcmc_runtime = pd.read_csv(require_file(mcmc_results_dir / "combined_runtime.csv"))
    runtime = abi_runtime[["dataset_id", "elapsed_sec"]].merge(
        mcmc_runtime[["dataset_id", "elapsed_sec"]],
        on="dataset_id",
        suffixes=("_abi", "_mcmc"),
        validate="one_to_one",
    )
    runtime["mcmc_to_abi_runtime_ratio"] = runtime["elapsed_sec_mcmc"] / runtime[
        "elapsed_sec_abi"
    ].replace(0.0, np.nan)
    runtime_summary = runtime.drop(columns="dataset_id").agg(["mean", "median"]).reset_index()
    diagnostics = pd.read_csv(require_file(mcmc_results_dir / "combined_chain_diagnostics.csv"))
    acceptance = pd.read_csv(require_file(mcmc_results_dir / "combined_acceptance.csv"))
    diagnostics = diagnostics[diagnostics["dataset_id"].isin(common_ids)].copy()
    acceptance = acceptance[acceptance["dataset_id"].isin(common_ids)].copy()
    decomposition_by_dataset, decomposition_summary = build_error_decomposition(
        parameter_pairs, diagnostics, abi_runtime
    )
    ess_sensitivity_rows: list[dict[str, float | int | str]] = []
    for ess_threshold in (0, 20, 50, 100):
        selected = decomposition_by_dataset.loc[
            decomposition_by_dataset["effective_sample_size"] >= ess_threshold
        ]
        for parameter in JOINT_PARAMETERS:
            group = selected.loc[selected["parameter"] == parameter]
            if group.empty:
                continue
            ess_sensitivity_rows.append(
                {
                    "ess_threshold": ess_threshold,
                    "parameter": parameter,
                    "n_parameter_datasets": int(len(group)),
                    "abi_total_mse_mc_corrected": float(
                        (group["abi_squared_error"] - group["abi_mc_variance"]).mean()
                    ),
                    "mcmc_target_mse_mc_corrected": float(
                        (group["mcmc_squared_error"] - group["mcmc_mc_variance"]).mean()
                    ),
                    "approximation_mse_mc_corrected": float(
                        (
                            group["approximation_squared_error"]
                            - group["abi_mc_variance"]
                            - group["mcmc_mc_variance"]
                        ).mean()
                    ),
                    "interaction_mc_corrected": float(
                        (group["interaction"] + 2.0 * group["mcmc_mc_variance"]).mean()
                    ),
                    "posterior_mean_correlation": float(
                        group["posterior_mean_abi"].corr(group["posterior_mean_mcmc"])
                    ),
                }
            )
    ess_sensitivity = pd.DataFrame(ess_sensitivity_rows)

    tables = {
        "parameter_method_summary": method_summary,
        "parameter_pairwise_by_dataset": parameter_pairs,
        "parameter_agreement_summary": parameter_agreement,
        "posterior_marginal_distances": marginal_distances,
        "posterior_marginal_distance_summary": marginal_summary,
        "posterior_joint_correlation_differences": joint_differences,
        "edge_probability_pairs": edge_pairs,
        "edge_agreement_by_dataset": edge_agreement,
        "edge_agreement_summary": edge_agreement_summary,
        "edge_metrics_by_dataset": edge_metrics,
        "edge_metric_method_summary": edge_metric_summary,
        "runtime_by_dataset": runtime,
        "runtime_summary": runtime_summary,
        "mcmc_chain_diagnostics": diagnostics,
        "mcmc_acceptance": acceptance,
        "error_decomposition_by_dataset": decomposition_by_dataset,
        "error_decomposition_summary": decomposition_summary,
        "error_decomposition_ess_sensitivity": ess_sensitivity,
    }
    for name, table in tables.items():
        table.to_csv(output_dir / f"{name}.csv", index=False)
    return tables


def _style_axis(axis, grid_axis: str = "y") -> None:
    axis.set_facecolor("white")
    axis.grid(axis=grid_axis, color="#d8d8d8", linewidth=0.7, alpha=0.65)
    axis.set_axisbelow(True)


def plot_parameter_summary(tables: dict[str, pd.DataFrame], output_dir: Path):
    summary = tables["parameter_method_summary"].copy()
    parameters = [
        parameter for parameter in JOINT_PARAMETERS if parameter in set(summary["parameter"])
    ]
    x = np.arange(len(parameters))
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), facecolor="white")
    for offset, (method, color) in enumerate(
        [("ABI", ABI_COLOR), ("Matched MCMC", MCMC_COLOR)]
    ):
        subset = summary.set_index(["method", "parameter"])
        values_mae = [subset.loc[(method, parameter), "mean_abs_error"] for parameter in parameters]
        values_coverage = [subset.loc[(method, parameter), "coverage_95"] for parameter in parameters]
        shift = (offset - 0.5) * width
        axes[0].bar(x + shift, values_mae, width=width, color=color, label=method)
        axes[1].bar(x + shift, values_coverage, width=width, color=color, label=method)
    axes[1].axhline(0.95, color=REFERENCE_COLOR, linestyle="--", linewidth=1.2)
    axes[0].set_ylabel("Mean absolute error")
    axes[1].set_ylabel("Empirical 95% coverage")
    for axis in axes:
        axis.set_xticks(x, [PARAMETER_LABELS[p] for p in parameters])
        _style_axis(axis)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    path = Path(output_dir) / "parameter_recovery_bars.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    fig.savefig(
        Path(output_dir) / "parameter_recovery_comparison.png",
        dpi=180,
        bbox_inches="tight",
    )
    return fig


def _joint_parameter_pairs(tables: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, list[str]]:
    pairs = tables["parameter_pairwise_by_dataset"]
    parameters = [
        parameter for parameter in JOINT_PARAMETERS if parameter in set(pairs["parameter"])
    ]
    return pairs, parameters


def _equalize_xy_limits(axis) -> None:
    lower = min(axis.get_xlim()[0], axis.get_ylim()[0])
    upper = max(axis.get_xlim()[1], axis.get_ylim()[1])
    padding = 0.03 * (upper - lower) if upper > lower else 0.1
    axis.set_xlim(lower - padding, upper + padding)
    axis.set_ylim(lower - padding, upper + padding)
    axis.plot(
        [lower - padding, upper + padding],
        [lower - padding, upper + padding],
        color=REFERENCE_COLOR,
        linestyle="--",
        linewidth=1.0,
        zorder=1,
    )


def plot_parameter_truth_scatter(tables: dict[str, pd.DataFrame], output_dir: Path):
    pairs, parameters = _joint_parameter_pairs(tables)
    fig, axes = plt.subplots(
        1, len(parameters), figsize=(4.4 * len(parameters), 4.8), facecolor="white"
    )
    axes = np.atleast_1d(axes)
    for index, (axis, parameter) in enumerate(zip(axes, parameters)):
        group = pairs.loc[pairs["parameter"] == parameter]
        truth = group["truth_abi"].to_numpy(dtype=float)
        for method, color, suffix in (
            ("ABI", ABI_COLOR, "abi"),
            ("Matched MCMC", MCMC_COLOR, "mcmc"),
        ):
            means = group[f"posterior_mean_{suffix}"].to_numpy(dtype=float)
            lower = group[f"lower_95_{suffix}"].to_numpy(dtype=float)
            upper = group[f"upper_95_{suffix}"].to_numpy(dtype=float)
            axis.errorbar(
                truth,
                means,
                yerr=np.vstack((np.maximum(means - lower, 0.0), np.maximum(upper - means, 0.0))),
                fmt="o",
                markersize=4.2,
                color=color,
                ecolor=color,
                elinewidth=0.65,
                alpha=0.48,
                label=method,
                zorder=2,
            )
        _equalize_xy_limits(axis)
        axis.set_title(PARAMETER_LABELS[parameter])
        axis.set_xlabel("Generating value")
        if index == 0:
            axis.set_ylabel("Posterior mean (95% interval)")
            axis.legend(frameon=False)
        _style_axis(axis, grid_axis="both")
    fig.tight_layout()
    path = Path(output_dir) / "parameter_recovery_truth_scatter.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_parameter_agreement_scatter(tables: dict[str, pd.DataFrame], output_dir: Path):
    pairs, parameters = _joint_parameter_pairs(tables)
    fig, axes = plt.subplots(
        1, len(parameters), figsize=(4.4 * len(parameters), 4.8), facecolor="white"
    )
    axes = np.atleast_1d(axes)
    for index, (axis, parameter) in enumerate(zip(axes, parameters)):
        group = pairs.loc[pairs["parameter"] == parameter]
        abi_mean = group["posterior_mean_abi"].to_numpy(dtype=float)
        mcmc_mean = group["posterior_mean_mcmc"].to_numpy(dtype=float)
        axis.errorbar(
            mcmc_mean,
            abi_mean,
            xerr=np.vstack(
                (
                    np.maximum(mcmc_mean - group["lower_95_mcmc"].to_numpy(dtype=float), 0.0),
                    np.maximum(group["upper_95_mcmc"].to_numpy(dtype=float) - mcmc_mean, 0.0),
                )
            ),
            yerr=np.vstack(
                (
                    np.maximum(abi_mean - group["lower_95_abi"].to_numpy(dtype=float), 0.0),
                    np.maximum(group["upper_95_abi"].to_numpy(dtype=float) - abi_mean, 0.0),
                )
            ),
            fmt="o",
            markersize=4.4,
            color=INTERACTION_COLOR,
            ecolor=INTERACTION_COLOR,
            elinewidth=0.65,
            alpha=0.45,
            zorder=2,
        )
        _equalize_xy_limits(axis)
        axis.set_title(PARAMETER_LABELS[parameter])
        axis.set_xlabel("Matched MCMC posterior mean")
        if index == 0:
            axis.set_ylabel("ABI posterior mean (95% intervals)")
        _style_axis(axis, grid_axis="both")
    fig.tight_layout()
    path = Path(output_dir) / "parameter_recovery_agreement_scatter.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def _add_grouped_boxplots(axis, pairs: pd.DataFrame, parameters: list[str], value: str) -> None:
    positions = np.arange(len(parameters), dtype=float)
    for offset, (method, color, suffix) in enumerate(
        (("ABI", ABI_COLOR, "abi"), ("Matched MCMC", MCMC_COLOR, "mcmc"))
    ):
        shift = (offset - 0.5) * 0.38
        data: list[np.ndarray] = []
        for parameter in parameters:
            group = pairs.loc[pairs["parameter"] == parameter]
            if value == "bias":
                values = group[f"posterior_mean_{suffix}"] - group[f"truth_{suffix}"]
            else:
                values = group[f"upper_95_{suffix}"] - group[f"lower_95_{suffix}"]
            data.append(values.dropna().to_numpy(dtype=float))
        boxplot = axis.boxplot(
            data,
            positions=positions + shift,
            widths=0.30,
            patch_artist=True,
            showfliers=False,
        )
        for box in boxplot["boxes"]:
            box.set(facecolor=color, alpha=0.72)
        for median in boxplot["medians"]:
            median.set(color=REFERENCE_COLOR)
    axis.set_xticks(positions, [PARAMETER_LABELS[p] for p in parameters])


def plot_parameter_bias_interval_boxplots(
    tables: dict[str, pd.DataFrame], output_dir: Path
):
    pairs, parameters = _joint_parameter_pairs(tables)
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.8), facecolor="white")
    _add_grouped_boxplots(axes[0], pairs, parameters, "bias")
    axes[0].axhline(0.0, color=REFERENCE_COLOR, linestyle="--", linewidth=1.0)
    axes[0].set_ylabel("Posterior-mean error")
    _add_grouped_boxplots(axes[1], pairs, parameters, "width")
    axes[1].set_ylabel("95% interval width")
    axes[1].legend(
        handles=[
            Patch(facecolor=ABI_COLOR, alpha=0.72, label="ABI"),
            Patch(facecolor=MCMC_COLOR, alpha=0.72, label="Matched MCMC"),
        ],
        frameon=False,
    )
    for axis in axes:
        _style_axis(axis)
    fig.tight_layout()
    path = Path(output_dir) / "parameter_bias_interval_boxplots.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_boundary_metric_bars(tables: dict[str, pd.DataFrame], output_dir: Path):
    summary = tables["edge_metric_method_summary"].set_index("method")
    methods = [method for method in ("ABI", "Matched MCMC") if method in summary.index]
    metrics = [
        ("auroc_mean", "AUROC"),
        ("average_precision_mean", "Average precision"),
        ("brier_mean", "Brier score"),
        ("sensitivity_mpm_mean", "MPM sensitivity"),
        ("specificity_mpm_mean", "MPM specificity"),
    ]
    colors = [ABI_COLOR if method == "ABI" else MCMC_COLOR for method in methods]
    fig, axes = plt.subplots(1, len(metrics), figsize=(17, 3.9), facecolor="white")
    for axis, (column, title) in zip(axes, metrics):
        values = [summary.loc[method, column] for method in methods]
        axis.bar(methods, values, color=colors, alpha=0.82)
        axis.set_title(title)
        axis.tick_params(axis="x", rotation=30)
        _style_axis(axis)
    axes[0].set_ylabel("Mean across datasets")
    fig.tight_layout()
    path = Path(output_dir) / "boundary_metric_bars.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_runtime_comparison(tables: dict[str, pd.DataFrame], output_dir: Path):
    runtime = tables["runtime_by_dataset"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), facecolor="white")
    axes[0].scatter(
        runtime["elapsed_sec_abi"],
        runtime["elapsed_sec_mcmc"],
        s=28,
        color=ABI_COLOR,
        alpha=0.70,
    )
    _equalize_xy_limits(axes[0])
    axes[0].set_xlabel("ABI elapsed seconds")
    axes[0].set_ylabel("Matched MCMC elapsed seconds")
    axes[1].hist(
        runtime["mcmc_to_abi_runtime_ratio"].dropna(),
        bins=15,
        color=ABI_COLOR,
        alpha=0.82,
        edgecolor="white",
    )
    axes[1].set_xlabel("Matched MCMC / ABI elapsed-time ratio")
    axes[1].set_ylabel("Dataset count")
    for axis in axes:
        _style_axis(axis)
    fig.tight_layout()
    path = Path(output_dir) / "runtime_comparison.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_error_decomposition(tables: dict[str, pd.DataFrame], output_dir: Path):
    summary = tables["error_decomposition_summary"].set_index("parameter")
    parameters = [parameter for parameter in JOINT_PARAMETERS if parameter in summary.index]
    x = np.arange(len(parameters))
    width = 0.19
    components = [
        ("abi_total_mse_mc_corrected", "ABI total MSE", REFERENCE_COLOR),
        ("mcmc_reference_mse_mc_corrected", "MCMC reference", MCMC_COLOR),
        ("approximation_mse_mc_corrected", "ABI-MCMC discrepancy", ABI_COLOR),
        ("interaction_mc_corrected", "Interaction", INTERACTION_COLOR),
    ]
    fig, ax = plt.subplots(figsize=(10.2, 4.8), facecolor="white")
    for index, (column, label, color) in enumerate(components):
        shift = (index - 1.5) * width
        ax.bar(
            x + shift,
            [summary.loc[parameter, column] for parameter in parameters],
            width=width,
            color=color,
            label=label,
        )
    ax.axhline(0.0, color=REFERENCE_COLOR, linewidth=0.9)
    ax.set_xticks(x, [PARAMETER_LABELS[parameter] for parameter in parameters])
    ax.set_ylabel("Mean squared error component")
    ax.legend(frameon=False, ncols=2)
    _style_axis(ax)
    fig.tight_layout()
    path = Path(output_dir) / "posterior_mean_error_decomposition.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_posterior_agreement(tables: dict[str, pd.DataFrame], output_dir: Path):
    distances = tables["posterior_marginal_distances"]
    pairs = tables["parameter_pairwise_by_dataset"]
    parameters = [
        parameter
        for parameter in JOINT_PARAMETERS
        if parameter in set(distances["parameter"])
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), facecolor="white")
    axes[0].boxplot(
        [distances.loc[distances["parameter"] == p, "wasserstein_over_mcmc_sd"].dropna() for p in parameters],
        tick_labels=[PARAMETER_LABELS[p] for p in parameters],
        patch_artist=True,
        boxprops={"facecolor": ABI_COLOR, "alpha": 0.72},
        medianprops={"color": REFERENCE_COLOR},
    )
    axes[0].set_ylabel("Wasserstein distance / MCMC posterior SD")
    axes[1].boxplot(
        [pairs.loc[pairs["parameter"] == p, "posterior_width_ratio_abi_to_mcmc"].dropna() for p in parameters],
        tick_labels=[PARAMETER_LABELS[p] for p in parameters],
        patch_artist=True,
        boxprops={"facecolor": MCMC_COLOR, "alpha": 0.72},
        medianprops={"color": REFERENCE_COLOR},
    )
    axes[1].axhline(1.0, color=REFERENCE_COLOR, linestyle="--", linewidth=1.2)
    axes[1].set_ylabel("ABI / MCMC 95% interval width")
    for axis in axes:
        _style_axis(axis)
    fig.tight_layout()
    path = Path(output_dir) / "posterior_distribution_agreement.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig


def plot_edge_agreement(tables: dict[str, pd.DataFrame], output_dir: Path):
    pairs = tables["edge_probability_pairs"]
    fig, ax = plt.subplots(figsize=(6.2, 5.8), facecolor="white")
    ax.scatter(
        pairs["boundary_prob_abi"],
        pairs["boundary_prob_mcmc"],
        s=8,
        alpha=0.12,
        color=ABI_COLOR,
        edgecolors="none",
    )
    ax.plot([0, 1], [0, 1], color=REFERENCE_COLOR, linestyle="--", linewidth=1.2)
    ax.set(xlabel="ABI boundary probability", ylabel="Matched MCMC boundary probability", xlim=(0, 1), ylim=(0, 1))
    _style_axis(ax, grid_axis="both")
    fig.tight_layout()
    path = Path(output_dir) / "edge_probability_agreement.png"
    fig.savefig(path, dpi=180, bbox_inches="tight")
    return fig
