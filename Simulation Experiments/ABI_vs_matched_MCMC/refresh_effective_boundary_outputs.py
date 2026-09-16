from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from effective_boundaries import effective_boundary_draws, effective_boundary_truth


DEFAULT_BANK = "benchmark_bank_seed123_n100"


def compute_auc(probability: np.ndarray, truth: np.ndarray) -> float:
    positive = np.flatnonzero(truth == 1)
    negative = np.flatnonzero(truth == 0)
    if positive.size == 0 or negative.size == 0:
        return float("nan")
    ranks = pd.Series(probability).rank(method="average").to_numpy()
    return float(
        (ranks[positive].sum() - positive.size * (positive.size + 1) / 2.0)
        / (positive.size * negative.size)
    )


def compute_average_precision(probability: np.ndarray, truth: np.ndarray) -> float:
    positive_total = int(np.sum(truth == 1))
    if positive_total == 0:
        return float("nan")
    order = np.argsort(-probability)
    ordered_truth = truth[order]
    true_positive = np.cumsum(ordered_truth == 1)
    precision = true_positive / np.arange(1, ordered_truth.size + 1)
    return float(np.sum(precision[ordered_truth == 1]) / positive_total)


def refresh_method(
    benchmark_dir: Path,
    manifest: pd.DataFrame,
    results_dir: Path,
    method: str,
    probability_column: str,
    median_column: str,
) -> pd.DataFrame:
    metric_rows: list[dict[str, float | int | str]] = []
    for row in manifest.itertuples(index=False):
        dataset_id = str(row.dataset_id)
        file_name = getattr(row, "file_name", f"{dataset_id}.npz")
        with np.load(benchmark_dir / "datasets" / file_name, allow_pickle=False) as archive:
            adjacency = np.asarray(archive["A"])
            dissimilarity = np.asarray(archive["Z"])
            filtered = np.asarray(archive["A_filtered"])
            edge_i = np.asarray(archive["edge_i"], dtype=int).reshape(-1)
            edge_j = np.asarray(archive["edge_j"], dtype=int).reshape(-1)
        dataset_output = results_dir / "per_dataset" / dataset_id
        draw_path = dataset_output / f"{dataset_id}_posterior_draws.csv.gz"
        eta_draws = pd.read_csv(draw_path, usecols=["eta"])["eta"].to_numpy()
        boundary_draws = effective_boundary_draws(
            eta_draws, adjacency, dissimilarity, edge_i, edge_j
        )
        probability = boundary_draws.mean(axis=0).astype(np.float32)
        median = (probability > 0.5).astype(np.int8)
        truth = effective_boundary_truth(filtered, edge_i, edge_j)
        count_draws = boundary_draws.sum(axis=1)
        count_interval = np.quantile(count_draws, [0.025, 0.975])
        edge_z = dissimilarity[edge_i, edge_j]

        edge_probabilities = pd.DataFrame(
            {
                "edge_index_1based": np.arange(edge_i.size) + 1,
                "node_i_1based": edge_i + 1,
                "node_j_1based": edge_j + 1,
                "edge_z": edge_z,
                "boundary_true": truth,
                "dataset_id": dataset_id,
                probability_column: probability,
                median_column: median,
            }
        )
        edge_probabilities.to_csv(
            dataset_output / f"{dataset_id}_edge_probabilities.csv", index=False
        )

        metrics = {
            "dataset_id": dataset_id,
            "edge_count": int(edge_i.size),
            "true_boundary_count": int(truth.sum()),
            "posterior_boundary_count_mpm": int(median.sum()),
            "auroc": compute_auc(probability, truth),
            "average_precision": compute_average_precision(probability, truth),
            "brier": float(np.mean((probability - truth) ** 2)),
            "sensitivity_mpm": float(median[truth == 1].mean())
            if np.any(truth == 1)
            else np.nan,
            "specificity_mpm": float((median[truth == 0] == 0).mean())
            if np.any(truth == 0)
            else np.nan,
            "boundary_count_mean_draws": float(count_draws.mean()),
            "boundary_count_lower_95": float(count_interval[0]),
            "boundary_count_upper_95": float(count_interval[1]),
            "boundary_count_truth_in_95": int(
                count_interval[0] <= truth.sum() <= count_interval[1]
            ),
        }
        pd.DataFrame([metrics]).to_csv(
            dataset_output / f"{dataset_id}_edge_metrics.csv", index=False
        )
        metric_rows.append(metrics)
        print(f"[{method}] refreshed {dataset_id}")

    combined = pd.DataFrame(metric_rows)
    combined.to_csv(results_dir / "combined_edge_metrics.csv", index=False)
    return combined


def refresh_r_inputs(benchmark_dir: Path, manifest: pd.DataFrame) -> pd.DataFrame:
    updated = manifest.copy()
    for row_index, row in enumerate(manifest.itertuples(index=False)):
        dataset_id = str(row.dataset_id)
        file_name = getattr(row, "file_name", f"{dataset_id}.npz")
        with np.load(benchmark_dir / "datasets" / file_name, allow_pickle=False) as archive:
            filtered = np.asarray(archive["A_filtered"])
            edge_i = np.asarray(archive["edge_i"], dtype=int).reshape(-1)
            edge_j = np.asarray(archive["edge_j"], dtype=int).reshape(-1)
        truth = effective_boundary_truth(filtered, edge_i, edge_j)
        edge_path = benchmark_dir / "r_inputs" / dataset_id / "edge_table.csv"
        edge_table = pd.read_csv(edge_path)
        edge_table["boundary_true"] = truth
        edge_table.to_csv(edge_path, index=False)
        updated.loc[row_index, "boundary_count"] = int(truth.sum())
    updated.to_csv(benchmark_dir / "benchmark_manifest.csv", index=False)
    return updated


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Refresh post-repair boundary outputs from saved posterior draws."
    )
    parser.add_argument("--benchmark-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    script_dir = Path(__file__).resolve().parent
    benchmark_dir = (
        args.benchmark_dir.resolve()
        if args.benchmark_dir is not None
        else script_dir / "datasets" / DEFAULT_BANK
    )
    manifest = pd.read_csv(benchmark_dir / "benchmark_manifest.csv")
    manifest = refresh_r_inputs(benchmark_dir, manifest)
    refresh_method(
        benchmark_dir,
        manifest,
        benchmark_dir / "abi_results_all100",
        "ABI",
        "boundary_prob_abi",
        "boundary_median_abi",
    )
    refresh_method(
        benchmark_dir,
        manifest,
        benchmark_dir / "matched_mcmc_results_all100",
        "Matched MCMC",
        "boundary_prob_mcmc",
        "boundary_median_mcmc",
    )


if __name__ == "__main__":
    main()
