from __future__ import annotations

import os
import sys
import time
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


os.environ.setdefault("KERAS_BACKEND", "tensorflow")
LOG2 = float(np.log(2.0))


def find_repository_root(start: Path | None = None) -> Path:
    current = (Path.cwd() if start is None else start).resolve()
    for candidate in (current, *current.parents):
        checkpoint = candidate / "Training" / "Checkpoints" / "poisson_dagar.keras"
        sampler = (
            candidate
            / "Simulation Experiments"
            / "ABI_vs_matched_MCMC"
            / "dagar_poisson_boundary_matched_mwg.cpp"
        )
        if checkpoint.exists() and sampler.exists():
            return candidate
    raise FileNotFoundError("Could not locate the ABI_poisson_regression repository root.")


def parse_dataset_ids(raw: object) -> list[str] | None:
    if raw is None:
        return None
    if isinstance(raw, (list, tuple)):
        parts = [str(piece).strip() for piece in raw]
    else:
        parts = [piece.strip() for piece in str(raw).split(",")]
    parts = [piece for piece in parts if piece]
    return parts or None


def summarize_draws(
    samples: np.ndarray, truth: float, dataset_id: str, parameter: str
) -> pd.DataFrame:
    q025, q50, q975 = np.quantile(samples, [0.025, 0.5, 0.975])
    posterior_mean = float(np.mean(samples))
    return pd.DataFrame(
        [
            {
                "dataset_id": dataset_id,
                "parameter": parameter,
                "posterior_mean": posterior_mean,
                "posterior_sd": float(np.std(samples, ddof=1)),
                "posterior_median": float(q50),
                "lower_95": float(q025),
                "upper_95": float(q975),
                "truth": float(truth),
                "bias_mean": posterior_mean - float(truth),
                "abs_error_mean": abs(posterior_mean - float(truth)),
                "covered_95": int(q025 <= truth <= q975),
            }
        ]
    )


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


def load_dataset(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def require_abi_dependencies():
    try:
        import bayesflow  # noqa: F401
        import keras
    except ImportError as exc:
        raise ImportError(
            "The ABI benchmark requires bayesflow, keras, and tensorflow.\n"
            f"Current interpreter: {sys.executable}\n"
            f"Python version: {sys.version.split()[0]}"
        ) from exc
    return keras


def set_global_seed(keras_module, seed: int) -> None:
    np.random.seed(seed)
    if hasattr(keras_module, "utils") and hasattr(keras_module.utils, "set_random_seed"):
        keras_module.utils.set_random_seed(seed)


def build_conditions(dataset: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    observations = np.asarray(dataset["obs"], dtype=np.float32)
    n = int(np.asarray(dataset["N"]).reshape(-1)[0])
    return {"obs": observations[np.newaxis, ...], "N": np.array([n], dtype=np.int32)}


def portable_path(path: Path, repository_root: Path) -> str:
    try:
        return path.resolve().relative_to(repository_root).as_posix()
    except ValueError:
        return str(path.resolve())


def save_config(path: Path, rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(list(rows)).to_csv(path, index=False)


def run_abi_benchmark(
    benchmark_dir: Path,
    results_dir: Path,
    checkpoint: Path,
    num_samples: int = 10_000,
    seed: int = 123,
    dataset_ids: object = None,
    max_datasets: int | None = None,
    save_draws: bool = True,
    quiet: bool = False,
) -> dict[str, object]:
    if num_samples <= 0:
        raise ValueError("num_samples must be positive.")
    if max_datasets is not None and max_datasets <= 0:
        raise ValueError("max_datasets must be positive when provided.")

    repository_root = find_repository_root()
    benchmark_dir = Path(benchmark_dir).resolve()
    results_dir = Path(results_dir).resolve()
    checkpoint = Path(checkpoint).resolve()
    manifest_path = benchmark_dir / "benchmark_manifest.csv"
    datasets_dir = benchmark_dir / "datasets"
    for required in (manifest_path, datasets_dir, checkpoint):
        if not required.exists():
            raise FileNotFoundError(f"Missing required input: {required}")

    results_dir.mkdir(parents=True, exist_ok=True)
    per_dataset_dir = results_dir / "per_dataset"
    per_dataset_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(manifest_path)
    selected_ids = parse_dataset_ids(dataset_ids)
    if selected_ids is not None:
        manifest = manifest.loc[manifest["dataset_id"].isin(selected_ids)].copy()
    if max_datasets is not None:
        manifest = manifest.head(max_datasets).copy()
    if manifest.empty:
        raise ValueError("No datasets selected.")

    keras_module = require_abi_dependencies()
    set_global_seed(keras_module, seed)
    if not quiet:
        print(f"Loading ABI checkpoint: {checkpoint}")
    workflow = keras_module.saving.load_model(str(checkpoint))

    save_config(
        results_dir / "abi_config.csv",
        [
            {
                "benchmark_dir": portable_path(benchmark_dir, repository_root),
                "results_dir": portable_path(results_dir, repository_root),
                "checkpoint": portable_path(checkpoint, repository_root),
                "mcmc_reference_sampler": "dagar_poisson_boundary_matched_mwg.cpp",
                "num_samples": int(num_samples),
                "seed": int(seed),
                "dataset_ids": "" if selected_ids is None else ",".join(selected_ids),
                "max_datasets": "" if max_datasets is None else int(max_datasets),
                "save_draws": bool(save_draws),
            }
        ],
    )

    parameter_rows: list[pd.DataFrame] = []
    edge_metric_rows: list[pd.DataFrame] = []
    runtime_rows: list[pd.DataFrame] = []

    for row_index, row in enumerate(manifest.itertuples(index=False), start=1):
        dataset_id = str(row.dataset_id)
        file_name = getattr(row, "file_name", f"{dataset_id}.npz")
        dataset = load_dataset(datasets_dir / file_name)
        conditions = build_conditions(dataset)
        m_value = float(np.asarray(dataset["M"]).reshape(-1)[0])
        edge_i = np.asarray(dataset["edge_i"], dtype=np.int32).reshape(-1)
        edge_j = np.asarray(dataset["edge_j"], dtype=np.int32).reshape(-1)
        edge_z = np.asarray(dataset["edge_z"], dtype=np.float32).reshape(-1)
        boundary_truth = np.asarray(dataset["edge_boundary_true"], dtype=np.int32).reshape(-1)

        if not quiet:
            print(
                f"[{row_index}/{len(manifest)}] Running ABI for {dataset_id} "
                f"(N={int(dataset['N'])}, edges={edge_i.size})"
            )
        start = time.perf_counter()
        posterior = workflow.sample(conditions=conditions, num_samples=num_samples)
        elapsed = time.perf_counter() - start

        beta0 = np.asarray(posterior["beta"], dtype=np.float32)[0, :, 0]
        sigma2_w = np.asarray(posterior["sigma2_w"], dtype=np.float32)[0, :, 0]
        eta_raw = np.asarray(posterior["eta_raw"], dtype=np.float32)[0, :, 0]
        rho = np.asarray(posterior["rho"], dtype=np.float32)[0, :, 0]
        eta = eta_raw * m_value
        draws = {
            "beta0": beta0,
            "sigma2_w": sigma2_w,
            "eta_raw": eta_raw,
            "eta": eta,
            "rho": rho,
        }
        truths = {
            "beta0": float(np.asarray(dataset["beta0_true"]).reshape(-1)[0]),
            "sigma2_w": float(np.asarray(dataset["sigma2_w_true"]).reshape(-1)[0]),
            "eta_raw": float(np.asarray(dataset["eta_raw_true"]).reshape(-1)[0]),
            "eta": float(np.asarray(dataset["eta_true"]).reshape(-1)[0]),
            "rho": float(np.asarray(dataset["rho_true"]).reshape(-1)[0]),
        }
        parameter_summary = pd.concat(
            [summarize_draws(draws[name], truths[name], dataset_id, name) for name in draws],
            ignore_index=True,
        )

        boundary_draws = eta[:, None] * edge_z[None, :] > LOG2
        boundary_probability = boundary_draws.mean(axis=0).astype(np.float32)
        boundary_median = (boundary_probability > 0.5).astype(np.int32)
        boundary_count_draws = boundary_draws.sum(axis=1)
        boundary_count_interval = np.quantile(boundary_count_draws, [0.025, 0.975])
        edge_probabilities = pd.DataFrame(
            {
                "edge_index_1based": np.arange(edge_i.size, dtype=np.int32) + 1,
                "node_i_1based": edge_i + 1,
                "node_j_1based": edge_j + 1,
                "edge_z": edge_z,
                "boundary_true": boundary_truth,
                "dataset_id": dataset_id,
                "boundary_prob_abi": boundary_probability,
                "boundary_median_abi": boundary_median,
            }
        )
        edge_metrics = pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "edge_count": int(edge_i.size),
                    "true_boundary_count": int(boundary_truth.sum()),
                    "posterior_boundary_count_mpm": int(boundary_median.sum()),
                    "auroc": compute_auc(boundary_probability, boundary_truth),
                    "average_precision": compute_average_precision(boundary_probability, boundary_truth),
                    "brier": float(np.mean((boundary_probability - boundary_truth) ** 2)),
                    "sensitivity_mpm": float(boundary_median[boundary_truth == 1].mean())
                    if np.any(boundary_truth == 1)
                    else np.nan,
                    "specificity_mpm": float((boundary_median[boundary_truth == 0] == 0).mean())
                    if np.any(boundary_truth == 0)
                    else np.nan,
                    "boundary_count_mean_draws": float(boundary_count_draws.mean()),
                    "boundary_count_lower_95": float(boundary_count_interval[0]),
                    "boundary_count_upper_95": float(boundary_count_interval[1]),
                    "boundary_count_truth_in_95": int(
                        boundary_count_interval[0]
                        <= boundary_truth.sum()
                        <= boundary_count_interval[1]
                    ),
                }
            ]
        )
        runtime = pd.DataFrame(
            [
                {
                    "dataset_id": dataset_id,
                    "elapsed_sec": elapsed,
                    "n_saved_draws": int(beta0.size),
                    "seconds_per_saved_draw": elapsed / max(1, beta0.size),
                    "seconds_per_1000_saved_draws": 1000.0 * elapsed / max(1, beta0.size),
                }
            ]
        )

        dataset_output = per_dataset_dir / dataset_id
        dataset_output.mkdir(parents=True, exist_ok=True)
        parameter_summary.to_csv(dataset_output / f"{dataset_id}_parameter_summary.csv", index=False)
        edge_probabilities.to_csv(dataset_output / f"{dataset_id}_edge_probabilities.csv", index=False)
        edge_metrics.to_csv(dataset_output / f"{dataset_id}_edge_metrics.csv", index=False)
        runtime.to_csv(dataset_output / f"{dataset_id}_runtime.csv", index=False)
        if save_draws:
            pd.DataFrame({"draw": np.arange(1, beta0.size + 1), **draws}).to_csv(
                dataset_output / f"{dataset_id}_posterior_draws.csv.gz",
                index=False,
                compression="gzip",
            )

        parameter_rows.append(parameter_summary)
        edge_metric_rows.append(edge_metrics)
        runtime_rows.append(runtime)

    parameter_df = pd.concat(parameter_rows, ignore_index=True)
    edge_metrics_df = pd.concat(edge_metric_rows, ignore_index=True)
    runtime_df = pd.concat(runtime_rows, ignore_index=True)
    parameter_df.to_csv(results_dir / "combined_parameter_summaries.csv", index=False)
    edge_metrics_df.to_csv(results_dir / "combined_edge_metrics.csv", index=False)
    runtime_df.to_csv(results_dir / "combined_runtime.csv", index=False)

    if not quiet:
        print(f"\nSaved ABI benchmark outputs to: {results_dir}")
    return {
        "parameter_df": parameter_df,
        "edge_metrics_df": edge_metrics_df,
        "runtime_df": runtime_df,
        "results_dir": results_dir,
    }
