from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path

import numpy as np


DEFAULT_BANK = "benchmark_bank_seed123_n100"


def find_repository_root(start: Path | None = None) -> Path:
    current = (Path.cwd() if start is None else start).resolve()
    for candidate in (current, *current.parents):
        if (candidate / "Training" / "Checkpoints" / "poisson_dagar.keras").exists():
            return candidate
    raise FileNotFoundError("Could not locate the ABI_poisson_regression repository root.")


def write_matrix(path: Path, matrix: np.ndarray) -> None:
    np.savetxt(path, matrix, delimiter=",", fmt="%.17g")


def write_r_bundle(path: Path, dataset: dict[str, np.ndarray]) -> None:
    path.mkdir(parents=True, exist_ok=True)

    y = np.asarray(dataset["y"]).reshape(-1)
    e = np.asarray(dataset["e"]).reshape(-1)
    x = np.asarray(dataset["x"]).reshape(-1)
    ordering = np.asarray(dataset["ordering"], dtype=np.int64).reshape(-1)

    with (path / "node_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_index_1based", "y", "e", "x", "ordering_1based"])
        for idx, (y_i, e_i, x_i, order_i) in enumerate(zip(y, e, x, ordering), start=1):
            writer.writerow([idx, int(y_i), float(e_i), float(x_i), int(order_i) + 1])

    edge_i = np.asarray(dataset["edge_i"], dtype=np.int64).reshape(-1)
    edge_j = np.asarray(dataset["edge_j"], dtype=np.int64).reshape(-1)
    edge_z = np.asarray(dataset["edge_z"], dtype=float).reshape(-1)
    filtered = np.asarray(dataset["A_filtered"])
    edge_truth = (filtered[edge_i, edge_j] < 0.5).astype(np.int64)
    with (path / "edge_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["edge_index_1based", "node_i_1based", "node_j_1based", "edge_z", "boundary_true"]
        )
        for idx, (i, j, z, truth) in enumerate(zip(edge_i, edge_j, edge_z, edge_truth), start=1):
            writer.writerow([idx, int(i) + 1, int(j) + 1, float(z), int(truth)])

    write_matrix(path / "A.csv", np.asarray(dataset["A"], dtype=float))
    write_matrix(path / "Z.csv", np.asarray(dataset["Z"], dtype=float))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare the fixed held-out dataset bank for the matched MCMC benchmark."
    )
    parser.add_argument("--bank-name", default=DEFAULT_BANK)
    parser.add_argument("--source-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repository_root = find_repository_root()
    experiment_dir = Path(__file__).resolve().parent
    output_dir = (
        args.output_dir.resolve()
        if args.output_dir is not None
        else experiment_dir / "datasets" / args.bank_name
    )
    source_dir = args.source_dir.resolve() if args.source_dir is not None else output_dir

    manifest_path = source_dir / "benchmark_manifest.csv"
    source_datasets_dir = source_dir / "datasets"
    if not manifest_path.exists() or not source_datasets_dir.exists():
        raise FileNotFoundError(
            "The source benchmark bank is incomplete. Expected the manifest and datasets under "
            f"{source_dir}"
        )

    output_datasets_dir = output_dir / "datasets"
    r_inputs_dir = output_dir / "r_inputs"
    output_datasets_dir.mkdir(parents=True, exist_ok=True)
    r_inputs_dir.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))

    for idx, row in enumerate(rows, start=1):
        dataset_id = row["dataset_id"]
        file_name = row.get("file_name") or f"{dataset_id}.npz"
        source_path = source_datasets_dir / file_name
        target_path = output_datasets_dir / file_name
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source dataset: {source_path}")

        if source_path.resolve() != target_path.resolve():
            shutil.copy2(source_path, target_path)
        with np.load(source_path, allow_pickle=False) as archive:
            dataset = {key: archive[key] for key in archive.files}
        write_r_bundle(r_inputs_dir / dataset_id, dataset)
        print(f"[{idx:>3d}/{len(rows)}] prepared {dataset_id}")

    output_manifest_path = output_dir / "benchmark_manifest.csv"
    if manifest_path.resolve() != output_manifest_path.resolve():
        shutil.copy2(manifest_path, output_manifest_path)
    source_config_path = source_dir / "benchmark_config.json"
    output_config_path = output_dir / "benchmark_config.json"
    if source_config_path.exists() and source_config_path.resolve() != output_config_path.resolve():
        shutil.copy2(source_config_path, output_config_path)
    try:
        source_bank = source_dir.relative_to(repository_root).as_posix()
    except ValueError:
        source_bank = str(source_dir)
    config = {
        "source_bank": source_bank,
        "bank_name": args.bank_name,
        "num_datasets": len(rows),
        "dataset_dir": "datasets",
        "r_input_dir": "r_inputs",
        "manifest_file": "benchmark_manifest.csv",
        "purpose": "prior- and construction-matched ABI versus MCMC comparison",
    }
    with (output_dir / "matched_benchmark_config.json").open("w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)
        handle.write("\n")

    print(f"\nPrepared matched benchmark bank: {output_dir}")


if __name__ == "__main__":
    main()
