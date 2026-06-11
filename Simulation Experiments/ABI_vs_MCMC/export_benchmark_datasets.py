from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np

from benchmark_utils import (
    N_MAX,
    N_MIN,
    ORDERING_MODE,
    dataset_summary,
    save_dataset_csv_bundle,
    save_dataset_npz,
    simulate_dataset,
    summary_to_dict,
    write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Export a benchmark bank of simulated datasets generated from the "
            "same ABI-DAGAR simulator used in the 'more general' notebook."
        )
    )
    parser.add_argument("--num-datasets", type=int, default=48, help="Number of datasets to export.")
    parser.add_argument("--seed", type=int, default=123, help="Master random seed.")
    parser.add_argument(
        "--run-name",
        type=str,
        default=None,
        help="Name of the exported bank. Defaults to benchmark_bank_seed<seed>_n<count>.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to ABI_vs_MCMC/datasets/<run-name>.",
    )
    parser.add_argument(
        "--ordering-mode",
        type=str,
        default=ORDERING_MODE,
        choices=["identity", "random"],
        help="Ordering rule used inside the simulator.",
    )
    parser.add_argument("--fixed-n", type=int, default=None, help="If provided, force all datasets to have this N.")
    parser.add_argument("--n-min", type=int, default=N_MIN, help="Minimum N when N is sampled.")
    parser.add_argument("--n-max", type=int, default=N_MAX, help="Maximum N when N is sampled.")
    return parser


def resolve_output_dir(script_dir: Path, run_name: str, requested: Path | None) -> Path:
    if requested is None:
        return script_dir / "datasets" / run_name
    if requested.is_absolute():
        return requested
    return (Path.cwd() / requested).resolve()


def write_manifest(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        raise ValueError("Cannot write an empty manifest.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.num_datasets <= 0:
        raise ValueError("--num-datasets must be positive.")
    if args.fixed_n is not None and args.fixed_n <= 0:
        raise ValueError("--fixed-n must be positive when provided.")
    if args.n_min <= 0 or args.n_max <= 0 or args.n_min > args.n_max:
        raise ValueError("--n-min and --n-max must be positive and satisfy n_min <= n_max.")

    run_name = args.run_name or f"benchmark_bank_seed{args.seed}_n{args.num_datasets}"
    script_dir = Path(__file__).resolve().parent
    output_dir = resolve_output_dir(script_dir, run_name, args.output_dir)
    datasets_dir = output_dir / "datasets"
    r_inputs_dir = output_dir / "r_inputs"
    datasets_dir.mkdir(parents=True, exist_ok=True)
    r_inputs_dir.mkdir(parents=True, exist_ok=True)

    master_rng = np.random.default_rng(args.seed)
    dataset_seeds = master_rng.integers(0, np.iinfo(np.uint32).max, size=args.num_datasets, dtype=np.uint32)

    manifest_rows: List[Dict[str, object]] = []

    for idx, dataset_seed in enumerate(dataset_seeds):
        dataset_id = f"dataset_{idx:04d}"
        file_name = f"{dataset_id}.npz"
        dataset_path = datasets_dir / file_name

        rng = np.random.default_rng(int(dataset_seed))
        dataset = simulate_dataset(
            rng=rng,
            n=args.fixed_n,
            ordering_mode=args.ordering_mode,
            n_min=args.n_min,
            n_max=args.n_max,
        )

        save_dataset_npz(dataset_path, dataset)

        summary = dataset_summary(
            dataset=dataset,
            dataset_id=dataset_id,
            file_name=file_name,
            dataset_seed=int(dataset_seed),
        )
        save_dataset_csv_bundle(r_inputs_dir / dataset_id, dataset, summary)
        manifest_rows.append(summary_to_dict(summary))

        print(
            f"[{idx + 1:>3d}/{args.num_datasets}] exported {dataset_id} "
            f"(N={summary.N}, edges={summary.edge_count}, boundaries={summary.boundary_count})"
        )

    write_manifest(output_dir / "benchmark_manifest.csv", manifest_rows)

    config = {
        "run_name": run_name,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "master_seed": int(args.seed),
        "num_datasets": int(args.num_datasets),
        "ordering_mode": args.ordering_mode,
        "fixed_n": args.fixed_n,
        "n_min": int(args.n_min),
        "n_max": int(args.n_max),
        "dataset_dir": str(datasets_dir),
        "r_input_dir": str(r_inputs_dir),
        "manifest_file": str(output_dir / "benchmark_manifest.csv"),
        "dataset_seeds": [int(x) for x in dataset_seeds.tolist()],
        "generator": "ABI-DAGAR simulator from 'poisson regression w spatial adj rndm graph transformer more general.ipynb'",
    }
    write_json(output_dir / "benchmark_config.json", config)

    print(f"\nSaved benchmark bank to: {output_dir}")
    print(f"Manifest: {output_dir / 'benchmark_manifest.csv'}")
    print(f"Config:   {output_dir / 'benchmark_config.json'}")


if __name__ == "__main__":
    main()
