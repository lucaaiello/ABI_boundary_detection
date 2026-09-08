"""Generate the empirical boundary-probability/residual-contrast figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
REAL_DATA_DIR = ROOT / "Real Data Analysis"
COMPARISON_DIR = REAL_DATA_DIR / "results_ABI_vs_CARBayes"
EXPLORATORY_DIR = REAL_DATA_DIR / "results_exploratory_analysis"
OUTPUT_PATH = (
    ROOT
    / "Figure Generation"
    / "Images"
    / "boundary_probability_residual_contrast.png"
)

APPLICATIONS = (
    ("glasgow", "Greater Glasgow", "#2a8192"),
    ("california", "California", "#d49b2a"),
    ("south_korea", "South Korea", "#bf6655"),
)


def canonical_id(value: object) -> str:
    """Normalize identifiers written by Python and R to the same key."""
    if pd.isna(value):
        return ""
    text = str(value).strip().lower()
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def load_edge_diagnostics(application: str) -> pd.DataFrame:
    edge_path = COMPARISON_DIR / application / "edge_comparison.csv"
    area_path = EXPLORATORY_DIR / f"{application}_harmonized_data.csv"

    edges = pd.read_csv(edge_path)
    areas = pd.read_csv(area_path)
    risk_by_id = {
        canonical_id(area_id): float(risk)
        for area_id, risk in zip(areas["area_id"], areas["crude_log_risk"])
    }

    from_risk = edges["from_id"].map(lambda value: risk_by_id.get(canonical_id(value)))
    to_risk = edges["to_id"].map(lambda value: risk_by_id.get(canonical_id(value)))
    diagnostics = pd.DataFrame(
        {
            "boundary_probability": pd.to_numeric(
                edges["boundary_prob_abi"], errors="coerce"
            ),
            "residual_contrast": np.abs(from_risk - to_risk),
        }
    ).dropna()

    if len(diagnostics) != len(edges):
        missing = len(edges) - len(diagnostics)
        raise ValueError(f"{application}: {missing} edges could not be matched to areas")
    return diagnostics


def binned_summary(data: pd.DataFrame) -> pd.DataFrame:
    bins = np.linspace(0.0, 1.0, 11)
    groups = pd.cut(
        data["boundary_probability"], bins=bins, include_lowest=True, duplicates="drop"
    )
    summary = (
        data.assign(probability_bin=groups)
        .groupby("probability_bin", observed=True)
        .agg(
            x=("boundary_probability", "mean"),
            y=("residual_contrast", "mean"),
            y_sd=("residual_contrast", "std"),
            count=("residual_contrast", "size"),
        )
        .reset_index(drop=True)
    )
    summary["y_sd"] = summary["y_sd"].fillna(0.0)
    return summary


def apply_axis_style(axis: plt.Axes) -> None:
    axis.set_facecolor("#fbfcfc")
    axis.grid(True, color="#8b969c", alpha=0.22, linewidth=0.8)
    axis.set_axisbelow(True)
    for spine in axis.spines.values():
        spine.set_visible(True)
        spine.set_color("#5f6b70")
        spine.set_linewidth(0.9)
    axis.tick_params(labelsize=10)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 11,
            "axes.titlesize": 13,
            "axes.titleweight": "bold",
        }
    )

    figure, axes = plt.subplots(1, 3, figsize=(12.35, 3.5), sharex=True)
    for index, (application, label, color) in enumerate(APPLICATIONS):
        axis = axes[index]
        data = load_edge_diagnostics(application)
        summary = binned_summary(data)
        correlation = spearmanr(
            data["boundary_probability"], data["residual_contrast"]
        ).statistic

        axis.scatter(
            data["boundary_probability"],
            data["residual_contrast"],
            s=14,
            color=color,
            alpha=0.27,
            edgecolors="none",
            rasterized=True,
        )
        axis.errorbar(
            summary["x"],
            summary["y"],
            yerr=summary["y_sd"],
            color="#26343b",
            linewidth=1.5,
            marker="o",
            markersize=4.8,
            markerfacecolor="#f7f9f9",
            markeredgecolor="#26343b",
            capsize=3,
            zorder=4,
        )
        axis.axvline(0.5, color="#6d777c", linestyle="--", linewidth=1.1, alpha=0.8)
        axis.text(
            0.97,
            0.955,
            rf"Spearman $r_s={correlation:.3f}$",
            transform=axis.transAxes,
            ha="right",
            va="top",
            fontsize=10.5,
            color="#46545b",
        )
        axis.set_title(label, pad=7)
        axis.set_xlim(-0.025, 1.025)
        axis.set_xlabel(r"ABI posterior boundary probability, $p_{ij}$")
        apply_axis_style(axis)

    axes[0].set_ylabel(r"Observed residual contrast, $|r_i-r_j|$")
    figure.subplots_adjust(left=0.063, right=0.994, bottom=0.18, top=0.89, wspace=0.28)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, dpi=450, facecolor="white")
    plt.close(figure)
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
