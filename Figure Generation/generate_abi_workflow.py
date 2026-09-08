"""Generate the conceptual ABI-DAGAR train, validate, and reuse workflow."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Polygon
import numpy as np
from scipy.spatial import Delaunay


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = (
    ROOT
    / "Figure Generation"
    / "Images"
    / "abi_train_validate_deploy_workflow.png"
)

INK = "#102a36"
TEAL = "#087582"
TEAL_LIGHT = "#e8f5f5"
AMBER = "#d68b00"
AMBER_LIGHT = "#fff6df"
CORAL = "#d62f1f"
CORAL_LIGHT = "#fff0ec"
EDGE = "#748b96"


def rounded_box(
    axis: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    edgecolor: str,
    facecolor: str = "white",
    linewidth: float = 1.4,
    radius: float = 0.012,
    linestyle: str = "-",
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
    )
    axis.add_patch(patch)
    return patch


def arrow(
    axis: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str,
    width: float = 1.8,
    connectionstyle: str = "arc3",
) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=width,
            color=color,
            connectionstyle=connectionstyle,
            shrinkA=0,
            shrinkB=0,
        )
    )


def stage_header(
    axis: plt.Axes, x: float, number: str, title: str, color: str, line_end: float
) -> None:
    axis.add_patch(Circle((x, 0.935), 0.019, facecolor=color, edgecolor="none"))
    axis.text(x, 0.935, number, ha="center", va="center", color="white", fontsize=13, weight="bold")
    axis.text(x + 0.028, 0.935, title, ha="left", va="center", color=color, fontsize=17, weight="bold")
    axis.plot([x + 0.17, line_end], [0.935, 0.935], color=color, linewidth=1.6)


def draw_graph(
    axis: plt.Axes,
    center: tuple[float, float],
    width: float,
    height: float,
    n_nodes: int,
    seed: int,
    accent: str,
    highlighted_edges: int = 2,
) -> None:
    rng = np.random.default_rng(seed)
    points = rng.uniform(0.05, 0.95, size=(n_nodes, 2))
    triangles = Delaunay(points).simplices
    edges: set[tuple[int, int]] = set()
    for triangle in triangles:
        for first, second in ((0, 1), (0, 2), (1, 2)):
            edge = tuple(sorted((int(triangle[first]), int(triangle[second]))))
            edges.add(edge)
    ordered_edges = sorted(edges)
    highlighted = set(rng.choice(len(ordered_edges), size=min(highlighted_edges, len(ordered_edges)), replace=False))

    x0, y0 = center
    coordinates = np.column_stack(
        [x0 + width * (points[:, 0] - 0.5), y0 + height * (points[:, 1] - 0.5)]
    )
    for edge_index, (first, second) in enumerate(ordered_edges):
        color = accent if edge_index in highlighted else EDGE
        linewidth = 1.8 if edge_index in highlighted else 0.65
        axis.plot(
            coordinates[[first, second], 0],
            coordinates[[first, second], 1],
            color=color,
            linewidth=linewidth,
            alpha=0.95,
            solid_capstyle="round",
            zorder=2,
        )
    axis.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        s=13,
        facecolor=INK,
        edgecolor="white",
        linewidth=0.45,
        zorder=3,
    )


def draw_network_block(axis: plt.Axes) -> None:
    rounded_box(axis, (0.025, 0.095), 0.218, 0.13, TEAL, facecolor="#f8fcfc")
    left_x = [0.043, 0.043, 0.043, 0.043]
    left_y = [0.13, 0.155, 0.18, 0.205]
    mid_x = [0.09, 0.09, 0.09]
    mid_y = [0.14, 0.165, 0.195]
    right_x = [0.13, 0.13, 0.13, 0.13]
    right_y = [0.13, 0.155, 0.18, 0.205]
    for x1, y1 in zip(left_x, left_y):
        for x2, y2 in zip(mid_x, mid_y):
            axis.plot([x1, x2], [y1, y2], color=EDGE, linewidth=0.45, alpha=0.75)
    for x1, y1 in zip(mid_x, mid_y):
        for x2, y2 in zip(right_x, right_y):
            axis.plot([x1, x2], [y1, y2], color=EDGE, linewidth=0.45, alpha=0.75)
    for x, y in list(zip(left_x, left_y)) + list(zip(mid_x, mid_y)) + list(zip(right_x, right_y)):
        axis.add_patch(Circle((x, y), 0.0045, facecolor="white", edgecolor=TEAL, linewidth=1.1))
    axis.text(0.151, 0.174, "SET TRANSFORMER", color=INK, fontsize=8.7, weight="bold", va="center")
    axis.text(0.151, 0.145, "+ CONDITIONAL FLOW", color=INK, fontsize=8.7, weight="bold", va="center")


def draw_posterior_operator(axis: plt.Axes) -> None:
    center = (0.35, 0.49)
    for radius, color, width in ((0.071, TEAL, 1.2), (0.058, "#51b9bf", 1.8), (0.044, TEAL, 1.0)):
        axis.add_patch(Circle(center, radius, facecolor="none", edgecolor=color, linewidth=width))
    axis.add_patch(Circle(center, 0.076, facecolor=TEAL_LIGHT, edgecolor=TEAL, linewidth=0.8, alpha=0.45))
    nodes = np.array(
        [
            [0.323, 0.515], [0.350, 0.545], [0.377, 0.512], [0.331, 0.475],
            [0.350, 0.490], [0.371, 0.472], [0.350, 0.445],
        ]
    )
    edges = ((0, 1), (1, 2), (0, 4), (2, 4), (3, 4), (4, 5), (3, 6), (5, 6))
    for first, second in edges:
        axis.plot(nodes[[first, second], 0], nodes[[first, second], 1], color=INK, linewidth=1.25)
    axis.scatter(nodes[:, 0], nodes[:, 1], s=23, facecolor="#d9f2f1", edgecolor=TEAL, linewidth=1.1, zorder=4)
    axis.text(center[0], 0.62, "ONE POSTERIOR", color=TEAL, ha="center", fontsize=11.5, weight="bold")
    axis.text(center[0], 0.59, "OPERATOR", color=TEAL, ha="center", fontsize=11.5, weight="bold")


def draw_validation(axis: plt.Axes) -> None:
    rounded_box(axis, (0.46, 0.16), 0.185, 0.65, AMBER, facecolor="#fffdfa", linestyle=(0, (3, 3)))
    axis.text(0.5525, 0.835, "FROZEN WEIGHTS", ha="center", color=AMBER, fontsize=10.8, weight="bold")

    shield = Polygon(
        [[0.515, 0.73], [0.5525, 0.77], [0.59, 0.73], [0.587, 0.62], [0.5525, 0.57], [0.518, 0.62]],
        closed=True,
        facecolor=AMBER_LIGHT,
        edgecolor=AMBER,
        linewidth=2.6,
    )
    axis.add_patch(shield)
    axis.add_patch(Arc((0.5525, 0.685), 0.033, 0.055, theta1=0, theta2=180, edgecolor=AMBER, linewidth=2.4))
    rounded_box(axis, (0.538, 0.637), 0.029, 0.043, AMBER, facecolor=AMBER, linewidth=1.0, radius=0.005)
    axis.add_patch(Circle((0.5525, 0.661), 0.0042, facecolor="white", edgecolor="none"))
    axis.plot([0.5525, 0.5525], [0.657, 0.647], color="white", linewidth=1.4)

    checklist = ("CALIBRATION", "BOUNDARY RECOVERY", "REPLICATION", "MCMC-DAGAR")
    for row, text in enumerate(checklist):
        y = 0.49 - row * 0.072
        rounded_box(axis, (0.475, y), 0.156, 0.058, "#a67a21", facecolor="white", linewidth=0.7, radius=0.004)
        rounded_box(axis, (0.484, y + 0.016), 0.012, 0.024, AMBER, facecolor=AMBER, linewidth=0.7, radius=0.002)
        axis.plot([0.487, 0.490, 0.494], [y + 0.028, y + 0.023, y + 0.034], color="white", linewidth=1.0)
        axis.text(0.505, y + 0.029, text, va="center", fontsize=7.8, color=INK, weight="bold")


def draw_gate(axis: plt.Axes) -> None:
    x0, y0 = 0.67, 0.49
    axis.plot([x0 - 0.026, x0 - 0.026], [y0 - 0.07, y0 + 0.07], color=INK, linewidth=3)
    axis.plot([x0 + 0.026, x0 + 0.026], [y0 - 0.07, y0 + 0.07], color=INK, linewidth=3)
    for offset in np.linspace(-0.021, 0.021, 7):
        axis.plot([x0 + offset, x0 + offset], [y0 - 0.058, y0 + 0.058], color=INK, linewidth=0.8)
    axis.add_patch(Circle((x0, y0), 0.021, facecolor=INK, edgecolor=AMBER, linewidth=1.2))
    axis.plot([x0 - 0.008, x0 - 0.001, x0 + 0.010], [y0, y0 - 0.008, y0 + 0.009], color=AMBER, linewidth=2.0)


def draw_reuse(axis: plt.Axes) -> None:
    axis.text(0.845, 0.865, "NO RETRAINING", ha="center", color=CORAL, fontsize=10.5, weight="bold")
    applications = (
        ("CALIFORNIA", 58, 0.735, 15, 21),
        ("GREATER GLASGOW", 134, 0.535, 18, 22),
        ("SOUTH KOREA", 241, 0.335, 22, 23),
    )
    for label, n_areas, y, nodes, seed in applications:
        rounded_box(axis, (0.735, y - 0.074), 0.235, 0.145, CORAL, facecolor="#fffdfc", linewidth=0.7)
        draw_graph(axis, (0.795, y), 0.105, 0.105, nodes, seed, CORAL, highlighted_edges=3)
        axis.text(0.855, y + 0.012, label, color=INK, fontsize=9.2, weight="bold", va="center")
        axis.text(0.855, y - 0.027, rf"$N={n_areas}$", color="#4d5d65", fontsize=9.2, va="center")

    rounded_box(axis, (0.75, 0.075), 0.205, 0.115, CORAL, facecolor=CORAL_LIGHT, linewidth=1.2)
    mini_nodes = np.array([[0.773, 0.13], [0.788, 0.148], [0.804, 0.13], [0.786, 0.105]])
    for first, second in ((0, 1), (1, 2), (0, 3), (2, 3)):
        color = CORAL if (first, second) == (1, 2) else INK
        axis.plot(mini_nodes[[first, second], 0], mini_nodes[[first, second], 1], color=color, linewidth=1.4)
    axis.scatter(mini_nodes[:, 0], mini_nodes[:, 1], s=24, facecolor=INK, edgecolor="white", linewidth=0.5, zorder=4)
    axis.text(0.875, 0.145, "JOINT POSTERIOR", ha="center", color=CORAL, fontsize=10.2, weight="bold")
    axis.text(0.875, 0.108, "EDGE PROBABILITIES", ha="center", color=INK, fontsize=9.4)


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    figure, axis = plt.subplots(figsize=(16, 6.56))
    figure.patch.set_facecolor("white")
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")

    rounded_box(axis, (0.012, 0.045), 0.412, 0.84, TEAL, facecolor=TEAL_LIGHT, linewidth=0.0, radius=0.018)
    rounded_box(axis, (0.44, 0.045), 0.225, 0.84, AMBER, facecolor=AMBER_LIGHT, linewidth=0.0, radius=0.018)
    rounded_box(axis, (0.69, 0.045), 0.298, 0.84, CORAL, facecolor=CORAL_LIGHT, linewidth=0.0, radius=0.018)
    for patch in axis.patches[:3]:
        patch.set_alpha(0.33)

    stage_header(axis, 0.032, "1", "TRAIN ONCE", TEAL, 0.415)
    stage_header(axis, 0.46, "2", "VALIDATE", AMBER, 0.65)
    stage_header(axis, 0.71, "3", "REUSE", CORAL, 0.978)

    axis.text(0.135, 0.84, "VARIABLE-SIZE GRAPHS", ha="center", color=INK, fontsize=10.5, weight="bold")
    axis.text(0.135, 0.803, r"$N=40,\ldots,300$", ha="center", color=INK, fontsize=10.5)
    graph_positions = ((0.06, 0.70), (0.135, 0.70), (0.21, 0.70), (0.06, 0.52), (0.135, 0.52), (0.21, 0.52))
    graph_sizes = (8, 10, 12, 13, 15, 17)
    for index, (position, n_nodes) in enumerate(zip(graph_positions, graph_sizes)):
        draw_graph(axis, position, 0.075, 0.105, n_nodes, 100 + index, CORAL, highlighted_edges=2)

    draw_network_block(axis)
    draw_posterior_operator(axis)
    arrow(axis, (0.243, 0.16), (0.285, 0.45), TEAL, connectionstyle="angle3,angleA=0,angleB=90")
    arrow(axis, (0.425, 0.49), (0.455, 0.49), AMBER)
    draw_validation(axis)
    arrow(axis, (0.647, 0.49), (0.655, 0.49), CORAL)
    draw_gate(axis)
    arrow(axis, (0.698, 0.49), (0.726, 0.49), CORAL)
    draw_reuse(axis)
    arrow(axis, (0.85, 0.255), (0.85, 0.195), CORAL)

    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PATH, dpi=300, facecolor="white")
    plt.close(figure)
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
