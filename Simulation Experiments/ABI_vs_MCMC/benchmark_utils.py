from __future__ import annotations

import json
import csv
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict

import numpy as np
from scipy.spatial import Delaunay


LOG2 = float(np.log(2.0))
ETA_RAW_LOWER = 0.0
ORDERING_MODE = "identity"
N_MIN = 40
N_MAX = 300
EXPOSURE_LOG_LOWER = float(np.log(2.0))
EXPOSURE_LOG_UPPER = float(np.log(30000.0))
POISSON_LAMBDA_MIN = 1e-2
POISSON_LAMBDA_MAX = 1e6


@dataclass
class DatasetSummary:
    dataset_id: str
    file_name: str
    dataset_seed: int
    N: int
    edge_count: int
    filtered_edge_count: int
    thresholded_edge_count: int
    boundary_count: int
    avg_degree: float
    avg_degree_filtered: float
    avg_degree_thresholded: float
    Z_median: float
    M: float
    beta0_true: float
    sigma2_w_true: float
    eta_raw_true: float
    eta_true: float
    rho_true: float


def random_adjacency(n: int, rng: np.random.Generator) -> np.ndarray:
    """Generate the same style of planar graph used in the notebook simulator."""
    points = rng.uniform(0.0, 10.0, size=(n, 2))
    tri = Delaunay(points)
    a = np.zeros((n, n), dtype=np.int32)
    for simplex in tri.simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                n1 = int(simplex[i])
                n2 = int(simplex[j])
                a[n1, n2] = 1
                a[n2, n1] = 1
    return a


def dagar_factors(a: np.ndarray, rho: float, ordering: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mirror the notebook DAGAR factor construction."""
    n = a.shape[0]
    rho2 = rho ** 2
    inv_order = np.argsort(ordering)

    b = np.zeros((n, n), dtype=np.float32)
    lam = np.zeros(n, dtype=np.float32)

    for pos in range(n):
        i = int(ordering[pos])
        preds = [int(ordering[q]) for q in range(pos) if a[i, int(ordering[q])] == 1]
        n_lt = len(preds)
        denom = 1.0 + max(n_lt - 1, 0) * rho2
        b_val = rho / denom if n_lt > 0 else 0.0
        for j in preds:
            b[pos, inv_order[j]] = b_val
        lam[pos] = denom / (1.0 - rho2)

    imb = np.eye(n, dtype=np.float32) - b
    return imb, lam


def get_ordering(a: np.ndarray, rng: np.random.Generator, ordering_mode: str = ORDERING_MODE) -> np.ndarray:
    n = a.shape[0]
    if ordering_mode == "random":
        return rng.permutation(n).astype(np.int32)
    if ordering_mode == "identity":
        return np.arange(n, dtype=np.int32)
    raise ValueError(f"Unknown ordering_mode: {ordering_mode}")


def repair_isolates_deterministic(a_filtered: np.ndarray, a: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Exact isolate repair rule from the notebook."""
    a_rep = a_filtered.copy().astype(np.float32)
    for i in range(a_rep.shape[0]):
        if a_rep[i].sum() == 0:
            neighbors = np.where(a[i] == 1)[0]
            if len(neighbors) > 0:
                j = neighbors[np.argmin(z[i, neighbors])]
                a_rep[i, j] = 1.0
                a_rep[j, i] = 1.0
    return a_rep


def masked_row_mean(values: np.ndarray, mask: np.ndarray) -> np.ndarray:
    mask_f = mask.astype(np.float32)
    denom = mask_f.sum(axis=1)
    denom_safe = np.where(denom == 0, 1.0, denom)
    out = (values * mask_f).sum(axis=1) / denom_safe
    out[denom == 0] = 0.0
    return out.astype(np.float32)


def safe_mean_1d(values: np.ndarray) -> float:
    return float(values.mean()) if values.size > 0 else 0.0


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float32)
    y = np.asarray(y, dtype=np.float32)
    x_c = x - x.mean()
    y_c = y - y.mean()
    denom = np.sqrt(np.mean(x_c ** 2) * np.mean(y_c ** 2))
    if denom < 1e-8:
        return 0.0
    return float(np.mean(x_c * y_c) / denom)


def observed_eta_signal_metrics(
    y: np.ndarray,
    e: np.ndarray,
    a: np.ndarray,
    z: np.ndarray,
    z_median: float,
    r_lag_all: np.ndarray | None = None,
) -> Dict[str, float]:
    a_bool = a == 1
    edge_i, edge_j = np.where(np.triu(a_bool, 1))

    zero_metrics = dict(
        edge_absdiff_low=0.0,
        edge_absdiff_mid=0.0,
        edge_absdiff_high=0.0,
        edge_concord_low=0.0,
        edge_concord_mid=0.0,
        edge_concord_high=0.0,
        edge_absdiff_slope=0.0,
        edge_absdiff_gap=0.0,
        edge_concord_gap=0.0,
        edge_corr_all=0.0,
        lag_corr_all=0.0,
        lag_slope_all=0.0,
        edge_semivar_all=0.0,
        local_moran_mean=0.0,
    )
    if edge_i.size == 0:
        return zero_metrics

    r = (np.log(y + 0.5) - np.log(e)).astype(np.float32)
    z_rel = (z / z_median).astype(np.float32)

    z_edge = z_rel[edge_i, edge_j]
    absdiff_edge = np.abs(r[edge_i] - r[edge_j]).astype(np.float32)
    sqdiff_edge = ((r[edge_i] - r[edge_j]) ** 2).astype(np.float32)

    r_centered = (r - r.mean()).astype(np.float32)
    var_r = float(np.mean(r_centered ** 2))
    var_r_safe = max(var_r, 1e-8)

    concord_edge = (r_centered[edge_i] * r_centered[edge_j]).astype(np.float32)

    low = z_edge <= 0.75
    mid = (z_edge > 0.75) & (z_edge <= 1.25)
    high = z_edge > 1.25

    z_bar = float(z_edge.mean())
    absdiff_bar = float(absdiff_edge.mean())
    var_z = float(np.mean((z_edge - z_bar) ** 2))

    if var_z < 1e-8:
        edge_absdiff_slope = 0.0
    else:
        edge_absdiff_slope = float(np.mean((z_edge - z_bar) * (absdiff_edge - absdiff_bar)) / var_z)

    if r_lag_all is None:
        degree = a_bool.sum(axis=1).astype(np.float32)
        degree_safe = np.where(degree == 0, 1.0, degree)
        w = a / degree_safe[:, None]
        r_lag_all = (w @ r).astype(np.float32)

    r_lag_centered = (r_lag_all - r_lag_all.mean()).astype(np.float32)
    lag_slope_all = float(np.mean(r_centered * r_lag_centered) / var_r_safe)
    lag_corr_all = safe_corr(r, r_lag_all)

    local_moran = (r_centered * r_lag_all / var_r_safe).astype(np.float32)
    edge_corr_all = float(np.mean(concord_edge) / var_r_safe)
    edge_semivar_all = float(0.5 * np.mean(sqdiff_edge))

    return dict(
        edge_absdiff_low=safe_mean_1d(absdiff_edge[low]),
        edge_absdiff_mid=safe_mean_1d(absdiff_edge[mid]),
        edge_absdiff_high=safe_mean_1d(absdiff_edge[high]),
        edge_concord_low=safe_mean_1d(concord_edge[low]),
        edge_concord_mid=safe_mean_1d(concord_edge[mid]),
        edge_concord_high=safe_mean_1d(concord_edge[high]),
        edge_absdiff_slope=edge_absdiff_slope,
        edge_absdiff_gap=safe_mean_1d(absdiff_edge[high]) - safe_mean_1d(absdiff_edge[low]),
        edge_concord_gap=safe_mean_1d(concord_edge[low]) - safe_mean_1d(concord_edge[high]),
        edge_corr_all=edge_corr_all,
        lag_corr_all=lag_corr_all,
        lag_slope_all=lag_slope_all,
        edge_semivar_all=edge_semivar_all,
        local_moran_mean=float(local_moran.mean()),
    )


def build_observed_features(
    x: np.ndarray,
    y: np.ndarray,
    e: np.ndarray,
    a: np.ndarray,
    z: np.ndarray,
    z_median: float,
    m: float,
) -> np.ndarray:
    a_bool = a == 1

    x = x.astype(np.float32)
    y = y.astype(np.float32)
    e = e.astype(np.float32)

    log_y = np.log1p(y).astype(np.float32)
    log_e = np.log(e).astype(np.float32)
    r = (np.log(y + 0.5) - log_e).astype(np.float32)

    r_centered = (r - r.mean()).astype(np.float32)
    var_r_safe = max(float(np.mean(r_centered ** 2)), 1e-8)

    degree = a_bool.sum(axis=1).astype(np.float32)
    degree_safe = np.where(degree == 0, 1.0, degree)

    z_rel = (z / z_median).astype(np.float32)
    neigh_r = np.broadcast_to(r[None, :], z_rel.shape).astype(np.float32)
    abs_r_diff = np.abs(r[:, None] - r[None, :]).astype(np.float32)

    low_mask = a_bool & (z_rel <= 0.75)
    mid_mask = a_bool & (z_rel > 0.75) & (z_rel <= 1.25)
    high_mask = a_bool & (z_rel > 1.25)

    r_lag_all = masked_row_mean(neigh_r, a_bool)
    absdiff_all = masked_row_mean(abs_r_diff, a_bool)

    r_lag_low = masked_row_mean(neigh_r, low_mask)
    r_lag_mid = masked_row_mean(neigh_r, mid_mask)
    r_lag_high = masked_row_mean(neigh_r, high_mask)

    absdiff_low = masked_row_mean(abs_r_diff, low_mask)
    absdiff_mid = masked_row_mean(abs_r_diff, mid_mask)
    absdiff_high = masked_row_mean(abs_r_diff, high_mask)

    prop_low = (low_mask.sum(axis=1) / degree_safe).astype(np.float32)
    prop_mid = (mid_mask.sum(axis=1) / degree_safe).astype(np.float32)
    prop_high = (high_mask.sum(axis=1) / degree_safe).astype(np.float32)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_z_rel = np.nanmean(np.where(a_bool, z_rel, np.nan), axis=1)
        max_z_rel = np.nanmax(np.where(a_bool, z_rel, np.nan), axis=1)

    mean_z_rel = np.nan_to_num(mean_z_rel, nan=0.0).astype(np.float32)
    max_z_rel = np.nan_to_num(max_z_rel, nan=0.0).astype(np.float32)

    edge_metrics = observed_eta_signal_metrics(
        y=y,
        e=e,
        a=a,
        z=z,
        z_median=z_median,
        r_lag_all=r_lag_all,
    )

    local_moran = (r_centered * r_lag_all / var_r_safe).astype(np.float32)
    local_semivar = (((r - r_lag_all) ** 2) / var_r_safe).astype(np.float32)

    obs = np.stack(
        [
            x,
            log_y,
            log_e,
            r,
            degree,
            r_lag_all,
            absdiff_all,
            r_lag_low,
            r_lag_mid,
            r_lag_high,
            absdiff_low,
            absdiff_mid,
            absdiff_high,
            prop_low,
            prop_mid,
            prop_high,
            mean_z_rel,
            max_z_rel,
            local_moran,
            local_semivar,
            np.full(len(x), m, dtype=np.float32),
            np.full(len(x), edge_metrics["edge_absdiff_slope"], dtype=np.float32),
            np.full(len(x), edge_metrics["edge_absdiff_gap"], dtype=np.float32),
            np.full(len(x), edge_metrics["edge_concord_gap"], dtype=np.float32),
            np.full(len(x), edge_metrics["edge_corr_all"], dtype=np.float32),
            np.full(len(x), edge_metrics["lag_corr_all"], dtype=np.float32),
            np.full(len(x), edge_metrics["lag_slope_all"], dtype=np.float32),
            np.full(len(x), edge_metrics["edge_semivar_all"], dtype=np.float32),
        ],
        axis=-1,
    )

    return obs.astype(np.float32)


def simulate_dataset(
    rng: np.random.Generator,
    n: int | None = None,
    ordering_mode: str = ORDERING_MODE,
    n_min: int = N_MIN,
    n_max: int = N_MAX,
) -> Dict[str, Any]:
    if n is None:
        n = int(rng.integers(n_min, n_max + 1))
    else:
        n = int(n)

    beta0 = float(rng.normal(0.0, 0.5))
    sigma2_w = float(abs(rng.normal(0.0, 0.5)))
    eta_raw = float(rng.uniform(ETA_RAW_LOWER, 1.0))
    rho = float(rng.uniform(0.0, 1.0))

    a = random_adjacency(n, rng).astype(np.float32)

    x = rng.normal(0.0, 1.0, size=n).astype(np.float32)
    z = np.abs(x[:, None] - x[None, :]).astype(np.float32)

    z_edges = z[a == 1]
    if z_edges.size == 0:
        z_median = 1.0
    else:
        z_median = float(np.median(z_edges) + 1e-8)

    m = float(LOG2 / z_median)
    eta = float(eta_raw * m)

    a_thresholded = a * ((z * eta) <= LOG2).astype(np.float32)
    a_filtered = repair_isolates_deterministic(a_thresholded, a, z).astype(np.float32)

    ordering = get_ordering(a, rng=rng, ordering_mode=ordering_mode)

    imb, dagar_lam = dagar_factors(a_filtered, rho, ordering)
    z_latent = rng.normal(size=n).astype(np.float32)
    rhs = (np.sqrt(sigma2_w) * z_latent / np.sqrt(dagar_lam)).astype(np.float32)
    w_true = np.linalg.solve(imb, rhs).astype(np.float32)
    w_true = (w_true - np.mean(w_true)).astype(np.float32)

    log_e = rng.uniform(EXPOSURE_LOG_LOWER, EXPOSURE_LOG_UPPER, size=n)
    e = np.exp(log_e).astype(np.float32)

    log_poisson_lam = np.log(e) + beta0 + w_true
    poisson_lam_true = np.clip(np.exp(log_poisson_lam), POISSON_LAMBDA_MIN, POISSON_LAMBDA_MAX).astype(np.float32)
    y = rng.poisson(poisson_lam_true).astype(np.int32)

    obs = build_observed_features(
        x=x,
        y=y.astype(np.float32),
        e=e,
        a=a,
        z=z,
        z_median=z_median,
        m=m,
    )

    edge_i, edge_j = np.where(np.triu(a == 1, 1))
    edge_boundary_true = ((z[edge_i, edge_j] * eta) > LOG2).astype(np.int8)

    return dict(
        N=np.int32(n),
        beta0_true=np.float32(beta0),
        sigma2_w_true=np.float32(sigma2_w),
        eta_raw_true=np.float32(eta_raw),
        eta_true=np.float32(eta),
        rho_true=np.float32(rho),
        Z_median=np.float32(z_median),
        M=np.float32(m),
        x=x.astype(np.float32),
        y=y.astype(np.int32),
        e=e.astype(np.float32),
        Z=z.astype(np.float32),
        A=a.astype(np.float32),
        A_thresholded=a_thresholded.astype(np.float32),
        A_filtered=a_filtered.astype(np.float32),
        ordering=ordering.astype(np.int32),
        dagar_lambda_true=dagar_lam.astype(np.float32),
        w_true=w_true.astype(np.float32),
        poisson_lam_true=poisson_lam_true.astype(np.float32),
        obs=obs.astype(np.float32),
        edge_i=edge_i.astype(np.int32),
        edge_j=edge_j.astype(np.int32),
        edge_z=z[edge_i, edge_j].astype(np.float32),
        edge_boundary_true=edge_boundary_true.astype(np.int8),
    )


def dataset_summary(
    dataset: Dict[str, Any],
    dataset_id: str,
    file_name: str,
    dataset_seed: int,
) -> DatasetSummary:
    a = dataset["A"]
    a_filtered = dataset["A_filtered"]
    a_thresholded = dataset["A_thresholded"]

    degree = a.sum(axis=1)
    degree_filtered = a_filtered.sum(axis=1)
    degree_thresholded = a_thresholded.sum(axis=1)

    return DatasetSummary(
        dataset_id=dataset_id,
        file_name=file_name,
        dataset_seed=int(dataset_seed),
        N=int(dataset["N"]),
        edge_count=int(dataset["edge_i"].size),
        filtered_edge_count=int(np.triu(a_filtered, 1).sum()),
        thresholded_edge_count=int(np.triu(a_thresholded, 1).sum()),
        boundary_count=int(dataset["edge_boundary_true"].sum()),
        avg_degree=float(degree.mean()),
        avg_degree_filtered=float(degree_filtered.mean()),
        avg_degree_thresholded=float(degree_thresholded.mean()),
        Z_median=float(dataset["Z_median"]),
        M=float(dataset["M"]),
        beta0_true=float(dataset["beta0_true"]),
        sigma2_w_true=float(dataset["sigma2_w_true"]),
        eta_raw_true=float(dataset["eta_raw_true"]),
        eta_true=float(dataset["eta_true"]),
        rho_true=float(dataset["rho_true"]),
    )


def save_dataset_npz(path: Path, dataset: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **dataset)


def save_dataset_csv_bundle(path: Path, dataset: Dict[str, Any], summary: DatasetSummary) -> None:
    """Save an R-friendly per-dataset bundle using only CSV files."""
    path.mkdir(parents=True, exist_ok=True)

    node_table_path = path / "node_table.csv"
    with node_table_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "node_index_1based",
                "x",
                "y",
                "e",
                "w_true",
                "poisson_lam_true",
            ]
        )
        for idx in range(int(dataset["N"])):
            writer.writerow(
                [
                    idx + 1,
                    float(dataset["x"][idx]),
                    int(dataset["y"][idx]),
                    float(dataset["e"][idx]),
                    float(dataset["w_true"][idx]),
                    float(dataset["poisson_lam_true"][idx]),
                ]
            )

    edge_table_path = path / "edge_table.csv"
    with edge_table_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "edge_index_1based",
                "node_i_1based",
                "node_j_1based",
                "edge_z",
                "boundary_true",
            ]
        )
        for idx in range(dataset["edge_i"].size):
            writer.writerow(
                [
                    idx + 1,
                    int(dataset["edge_i"][idx]) + 1,
                    int(dataset["edge_j"][idx]) + 1,
                    float(dataset["edge_z"][idx]),
                    int(dataset["edge_boundary_true"][idx]),
                ]
            )

    np.savetxt(path / "A.csv", dataset["A"], delimiter=",", fmt="%.0f")
    np.savetxt(path / "Z.csv", dataset["Z"], delimiter=",", fmt="%.8f")
    np.savetxt(path / "A_thresholded.csv", dataset["A_thresholded"], delimiter=",", fmt="%.0f")
    np.savetxt(path / "A_filtered.csv", dataset["A_filtered"], delimiter=",", fmt="%.0f")
    np.savetxt(path / "obs.csv", dataset["obs"], delimiter=",", fmt="%.8f")

    metadata_path = path / "metadata.csv"
    metadata = summary_to_dict(summary)
    with metadata_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(metadata.keys()))
        writer.writeheader()
        writer.writerow(metadata)


def write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)


def summary_to_dict(summary: DatasetSummary) -> Dict[str, Any]:
    return asdict(summary)
