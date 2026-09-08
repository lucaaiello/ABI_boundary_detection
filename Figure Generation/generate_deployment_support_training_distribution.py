"""Compare empirical inputs with the original ABI-DAGAR training distribution."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.spatial import Delaunay


ROOT = Path(__file__).resolve().parents[1]
REAL_DATA_DIR = ROOT / "Real Data Analysis"
DATA_DIR = REAL_DATA_DIR / "Data"
EXPLORATORY_DIR = REAL_DATA_DIR / "results_exploratory_analysis"
OUTPUT_PATH = (
    ROOT
    / "Figure Generation"
    / "Images"
    / "deployment_support_training_distribution.png"
)
METRICS_PATH = EXPLORATORY_DIR / "deployment_support_metrics_seed20260907_n3000.csv"

SEED = 20260907
N_CONFIGURATIONS = 3000
LOG2 = float(np.log(2.0))

APPLICATIONS = (
    ("california", "California", "#d49b2a", DATA_DIR / "adjacency_matrix_california.csv"),
    ("glasgow", "Glasgow", "#1f788b", DATA_DIR / "adjacency_matrix_glasgow.csv"),
    (
        "south_korea",
        "South Korea",
        "#bd5b49",
        DATA_DIR / "South_Korea" / "adjacency_matrix_south_korea.csv",
    ),
)

METRICS = (
    ("n_areas", r"Number of areas, $N$"),
    ("n_edges", r"Adjacent edges, $|\mathcal{E}|$"),
    ("mean_degree", "Mean degree"),
    ("median_degree", "Median degree"),
    ("degree_q90", "90th percentile of degree"),
    ("median_edge_dissimilarity", r"Median neighboring $|\Delta x|$"),
    ("scale_m", r"Graph-specific scale, $M$"),
    ("covariate_moran_i", r"Covariate Moran's $I$"),
    ("median_log_exposure", r"Median $\log(e_i)$"),
)


def adjacency_from_points(points: np.ndarray) -> np.ndarray:
    triangulation = Delaunay(points)
    adjacency = np.zeros((len(points), len(points)), dtype=np.int8)
    for simplex in triangulation.simplices:
        adjacency[simplex[0], simplex[1]] = 1
        adjacency[simplex[1], simplex[0]] = 1
        adjacency[simplex[0], simplex[2]] = 1
        adjacency[simplex[2], simplex[0]] = 1
        adjacency[simplex[1], simplex[2]] = 1
        adjacency[simplex[2], simplex[1]] = 1
    return adjacency


def summarize_inputs(
    adjacency: np.ndarray, covariate: np.ndarray, exposure: np.ndarray
) -> dict[str, float]:
    n_areas = len(covariate)
    edge_i, edge_j = np.where(np.triu(adjacency, 1) == 1)
    degree = adjacency.sum(axis=1).astype(float)

    covariate = np.asarray(covariate, dtype=float)
    covariate_z = (covariate - covariate.mean()) / covariate.std(ddof=0)
    edge_dissimilarity = np.abs(covariate_z[edge_i] - covariate_z[edge_j])
    median_dissimilarity = float(np.median(edge_dissimilarity))

    denominator = float(covariate_z @ covariate_z)
    moran_i = (
        n_areas
        / float(adjacency.sum())
        * float(covariate_z @ adjacency @ covariate_z)
        / denominator
    )
    return {
        "n_areas": float(n_areas),
        "n_edges": float(len(edge_i)),
        "mean_degree": float(degree.mean()),
        "median_degree": float(np.median(degree)),
        "degree_q90": float(np.quantile(degree, 0.90)),
        "median_edge_dissimilarity": median_dissimilarity,
        "scale_m": LOG2 / median_dissimilarity,
        "covariate_moran_i": moran_i,
        "median_log_exposure": float(np.median(np.log(exposure))),
    }


def simulate_training_inputs() -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for index in range(N_CONFIGURATIONS):
        n_areas = int(rng.integers(40, 301))
        points = rng.uniform(0.0, 10.0, size=(n_areas, 2))
        adjacency = adjacency_from_points(points)
        covariate = rng.normal(0.0, 1.0, size=n_areas)
        log_exposure = rng.uniform(np.log(2.0), np.log(30000.0), size=n_areas)
        row = summarize_inputs(adjacency, covariate, np.exp(log_exposure))
        row.update(source="training", application="", configuration=index)
        rows.append(row)
    return pd.DataFrame(rows)


def load_empirical_inputs() -> pd.DataFrame:
    rows = []
    for application, label, _, adjacency_path in APPLICATIONS:
        areas = pd.read_csv(EXPLORATORY_DIR / f"{application}_harmonized_data.csv")
        adjacency = pd.read_csv(adjacency_path).to_numpy(dtype=np.int8)
        if adjacency.shape != (len(areas), len(areas)):
            raise ValueError(f"{application}: area and adjacency dimensions differ")
        row = summarize_inputs(
            adjacency,
            areas["covariate"].to_numpy(dtype=float),
            areas["expected"].to_numpy(dtype=float),
        )
        row.update(source="empirical", application=label, configuration=np.nan)
        rows.append(row)
    return pd.DataFrame(rows)


def apply_axis_style(axis: plt.Axes) -> None:
    axis.set_facecolor("#fbfcfc")
    axis.grid(True, color="#889399", alpha=0.22, linewidth=0.8)
    axis.set_axisbelow(True)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("#5f6b70")
        spine.set_linewidth(0.9)
    axis.tick_params(labelsize=9.5)


def main() -> None:
    training = simulate_training_inputs()
    empirical = load_empirical_inputs()
    combined = pd.concat([training, empirical], ignore_index=True)
    METRICS_PATH.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(METRICS_PATH, index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
        }
    )
    figure, axes = plt.subplots(3, 3, figsize=(11.55, 8.25))
    for axis, (metric, title) in zip(axes.flat, METRICS):
        values = training[metric].to_numpy(dtype=float)
        q025, q975 = np.quantile(values, [0.025, 0.975])
        axis.axvspan(q025, q975, color="#e5efed", alpha=0.95, zorder=0)
        axis.hist(
            values,
            bins=30,
            density=True,
            color="#d5dfdd",
            edgecolor="#edf2f1",
            linewidth=0.35,
            alpha=0.95,
            zorder=1,
        )
        for application, label, color, _ in APPLICATIONS:
            value = float(empirical.loc[empirical["application"] == label, metric].iloc[0])
            axis.axvline(value, color=color, linewidth=1.6, zorder=3)
            axis.plot(
                value,
                0.955,
                marker="v",
                markersize=6,
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.55,
                transform=axis.get_xaxis_transform(),
                clip_on=False,
                zorder=4,
            )
        axis.set_title(title, pad=7)
        axis.set_ylabel("Training density", fontsize=9.5, color="#637078")
        apply_axis_style(axis)

    legend_handles = [
        Line2D([0], [0], color=color, linewidth=2.2, label=label)
        for _, label, color, _ in APPLICATIONS
    ]
    figure.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.01),
        fontsize=10,
    )
    figure.subplots_adjust(
        left=0.065, right=0.99, bottom=0.09, top=0.965, hspace=0.38, wspace=0.28
    )
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, dpi=450, facecolor="white")
    plt.close(figure)
    print(f"Saved metrics: {METRICS_PATH}")
    print(f"Saved figure:  {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
