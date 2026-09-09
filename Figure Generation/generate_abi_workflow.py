"""Generate the polished ABI-DAGAR train, validate, and reuse workflow.

The layout mirrors the manuscript asset while remaining deterministic and
reproducible. Application thumbnails use the actual areal graphs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import Delaunay


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_PATH = ROOT / "Figure Generation" / "Images" / "abi_train_validate_deploy_workflow.png"
DATA_DIR = ROOT / "Real Data Analysis" / "Data"

CANVAS_WIDTH = 1958
CANVAS_HEIGHT = 803
OUTPUT_DPI = 100

INK = "#082734"
TEAL = "#006f79"
TEAL_MID = "#24a6ad"
TEAL_LIGHT = "#bcebee"
AMBER = "#d78900"
AMBER_LIGHT = "#fff1c8"
CORAL = "#d92b18"
EDGE = "#76909a"
WHITE = "#ffffff"
PAPER = WHITE


@dataclass(frozen=True)
class ApplicationGraph:
    label: str
    n_areas: int
    gpkg_path: Path
    adjacency_path: Path
    covariate: str
    center: tuple[float, float]
    box: tuple[float, float]
    node_size: float


APPLICATIONS = (
    ApplicationGraph(
        "CALIFORNIA",
        58,
        DATA_DIR / "respiratory_data_california.gpkg",
        DATA_DIR / "adjacency_matrix_california.csv",
        "smoking",
        (1572, 585),
        (166, 142),
        11,
    ),
    ApplicationGraph(
        "GLASGOW",
        134,
        DATA_DIR / "respiratory_data_glasgow.gpkg",
        DATA_DIR / "adjacency_matrix_glasgow.csv",
        "incomedep",
        (1572, 427),
        (188, 116),
        7,
    ),
    ApplicationGraph(
        "SOUTH KOREA",
        241,
        DATA_DIR / "South_Korea" / "mortality_data_south_korea.gpkg",
        DATA_DIR / "South_Korea" / "adjacency_matrix_south_korea.csv",
        "smoking_pct",
        (1572, 267),
        (170, 154),
        5,
    ),
)


def blend(color: str, toward: str = WHITE, amount: float = 0.25):
    source = np.asarray(matplotlib.colors.to_rgb(color))
    target = np.asarray(matplotlib.colors.to_rgb(toward))
    return tuple((1.0 - amount) * source + amount * target)


def add_shadow(artist, offset=(2.2, -2.2), alpha=0.16):
    artist.set_path_effects(
        [
            path_effects.SimplePatchShadow(offset=offset, alpha=alpha, rho=0.98),
            path_effects.Normal(),
        ]
    )
    return artist


def rounded_box(
    axis,
    xy,
    width,
    height,
    edgecolor,
    facecolor=WHITE,
    linewidth=1.5,
    radius=12,
    linestyle="-",
    shadow=False,
    zorder=2,
):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle=f"round,pad=0.0,rounding_size={radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        linestyle=linestyle,
        zorder=zorder,
    )
    axis.add_patch(patch)
    if shadow:
        add_shadow(patch)
    return patch


def arrow(
    axis,
    start,
    end,
    color,
    linewidth=3.0,
    mutation_scale=18,
    connectionstyle="arc3",
    zorder=8,
):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation_scale,
        linewidth=linewidth,
        color=color,
        connectionstyle=connectionstyle,
        shrinkA=0,
        shrinkB=0,
        joinstyle="miter",
        zorder=zorder,
    )
    axis.add_patch(patch)
    return patch


def draw_badge(axis, center, number, color):
    x, y = center
    axis.add_patch(Circle((x + 3, y - 4), 27, facecolor="#000000", edgecolor="none", alpha=0.12, zorder=3))
    for radius, amount in ((27, 0.0), (23, 0.08), (18, 0.16)):
        axis.add_patch(Circle((x, y), radius, facecolor=blend(color, amount=amount), edgecolor="none", zorder=4))
    axis.add_patch(Arc((x - 3, y + 5), 38, 29, theta1=35, theta2=150, color=WHITE, alpha=0.22, linewidth=2.0, zorder=5))
    axis.text(x, y - 1, number, ha="center", va="center", color=WHITE, fontsize=24, weight="bold", zorder=6)


def draw_stage_headers(axis):
    specs = (
        ((65, 746), "1", "TRAIN ONCE", TEAL, 112, 846, True),
        ((924, 746), "2", "VALIDATE", AMBER, 971, 1457, True),
        ((1515, 746), "3", "REUSE", CORAL, 1562, 1904, False),
    )
    widths = {"TRAIN ONCE": 275, "VALIDATE": 214, "REUSE": 154}
    for center, number, title, color, title_x, line_end, has_arrow in specs:
        draw_badge(axis, center, number, color)
        axis.text(title_x, 747, title, ha="left", va="center", color=color, fontsize=28, weight="bold")
        line_start = title_x + widths[title]
        if has_arrow:
            arrow(axis, (line_start, 746), (line_end, 746), color, linewidth=2.7, mutation_scale=16)
        else:
            axis.plot([line_start, line_end], [746, 746], color=color, linewidth=2.7, zorder=4)


def simulated_graph(seed, n_nodes):
    rng = np.random.default_rng(seed)
    points = rng.uniform((0.04, 0.05), (0.96, 0.95), size=(n_nodes, 2))
    triangles = Delaunay(points).simplices
    edges = set()
    for triangle in triangles:
        for first, second in ((0, 1), (0, 2), (1, 2)):
            edges.add(tuple(sorted((int(triangle[first]), int(triangle[second])))))
    edges = sorted(edges)
    highlighted = set(rng.choice(len(edges), size=min(2 + n_nodes // 8, len(edges)), replace=False))
    return points, edges, highlighted


def draw_simulated_graph(axis, center, width, height, n_nodes, seed):
    points, edges, highlighted = simulated_graph(seed, n_nodes)
    coordinates = np.column_stack(
        [center[0] + width * (points[:, 0] - 0.5), center[1] + height * (points[:, 1] - 0.5)]
    )
    for index, (first, second) in enumerate(edges):
        is_boundary = index in highlighted
        axis.plot(
            coordinates[[first, second], 0],
            coordinates[[first, second], 1],
            color=CORAL if is_boundary else EDGE,
            linewidth=2.4 if is_boundary else 1.0,
            alpha=0.97,
            solid_capstyle="round",
            zorder=3 if is_boundary else 2,
        )
    axis.scatter(coordinates[:, 0], coordinates[:, 1], s=23, facecolor=TEAL, edgecolor=INK, linewidth=0.6, zorder=4)


def draw_training_graphs(axis):
    axis.text(254, 676, "VARIABLE-SIZE GRAPHS", ha="center", va="center", color=INK, fontsize=18, weight="bold")
    axis.text(254, 641, r"$N=40,\ldots,300$", ha="center", va="center", color=INK, fontsize=17)
    positions = (
        (88, 570), (195, 570), (303, 570), (410, 570),
        (95, 420), (250, 420), (405, 420),
        (112, 270), (315, 270),
    )
    node_counts = (8, 9, 10, 11, 12, 13, 14, 15, 17)
    widths = (84, 91, 93, 94, 112, 118, 112, 145, 142)
    for index, (center, n_nodes, width) in enumerate(zip(positions, node_counts, widths)):
        draw_simulated_graph(axis, center, width, 75, n_nodes, 1800 + index)
    for index in range(6):
        axis.add_patch(Circle((210 + index * 22, 190), 4.1, facecolor="#d5dcdf" if index == 0 else "#a6b2b8", edgecolor="none", zorder=3))


def draw_network_block(axis):
    rounded_box(axis, (34, 42), 441, 105, TEAL, linewidth=2.0, radius=8, shadow=True)
    axis.text(118, 133, "SET TRANSFORMER", color=INK, fontsize=10.5, weight="bold", ha="center", va="center")
    axis.text(354, 133, "NORMALIZING FLOW", color=INK, fontsize=10.5, weight="bold", ha="center", va="center")

    set_nodes = np.array([[59, 72], [59, 96], [59, 120]])
    attention_nodes = np.array([[103, 80], [103, 112]])
    summary_node = np.array([[155, 96]])
    for source in set_nodes:
        for target in attention_nodes:
            axis.plot([source[0], target[0]], [source[1], target[1]], color=EDGE, linewidth=0.75, alpha=0.72, zorder=3)
    for source in attention_nodes:
        for target in summary_node:
            axis.plot([source[0], target[0]], [source[1], target[1]], color=EDGE, linewidth=0.75, alpha=0.72, zorder=3)
    for points in (set_nodes, attention_nodes):
        axis.scatter(points[:, 0], points[:, 1], s=53, facecolor=WHITE, edgecolor=TEAL, linewidth=1.8, zorder=5)
    axis.scatter(summary_node[:, 0], summary_node[:, 1], s=72, facecolor=TEAL_LIGHT, edgecolor=TEAL, linewidth=2.0, zorder=6)

    context_node = (245, 111)
    arrow(axis, (163, 98), (238, 110), TEAL, linewidth=1.7, mutation_scale=12)
    axis.text(199, 114, "SUMMARY", color=TEAL, fontsize=7.2, weight="bold", ha="center", va="bottom")
    axis.scatter([context_node[0]], [context_node[1]], s=60, facecolor=TEAL_LIGHT, edgecolor=TEAL, linewidth=1.8, zorder=6)

    base_node = (245, 80)
    axis.scatter([base_node[0]], [base_node[1]], s=49, facecolor=WHITE, edgecolor=INK, linewidth=1.5, zorder=6)
    axis.text(base_node[0], 55, r"$z_0$", color=INK, fontsize=9.0, ha="center", va="center")
    flow_centers = (286, 326, 366)
    previous_x = base_node[0] + 5
    for index, x in enumerate(flow_centers, start=1):
        arrow(axis, (previous_x, 80), (x - 13, 80), TEAL, linewidth=1.4, mutation_scale=9, zorder=5)
        rounded_box(axis, (x - 12, 66), 24, 28, TEAL, facecolor=blend(TEAL_LIGHT, amount=0.35), linewidth=1.2, radius=4, zorder=6)
        axis.text(x, 80, rf"$f_{index}$", color=INK, fontsize=8.0, ha="center", va="center", zorder=7)
        axis.plot([x, x], [111, 95], color=TEAL, linewidth=1.1, linestyle=(0, (2, 2)), zorder=4)
        previous_x = x + 13
    axis.plot([context_node[0], flow_centers[-1]], [111, 111], color=TEAL, linewidth=1.2, zorder=4)

    posterior_node = (423, 80)
    arrow(axis, (previous_x, 80), (posterior_node[0] - 7, 80), TEAL, linewidth=1.4, mutation_scale=9, zorder=5)
    axis.scatter([posterior_node[0]], [posterior_node[1]], s=66, facecolor=INK, edgecolor=TEAL, linewidth=1.4, zorder=6)
    axis.text(posterior_node[0], 55, r"$\theta$", color=INK, fontsize=9.0, ha="center", va="center")
    axis.plot([423, 475, 511, 511], [80, 80, 80, 353], color=TEAL, linewidth=3.0, solid_capstyle="round", zorder=7)
    arrow(axis, (511, 353), (598, 353), TEAL, linewidth=3.0, mutation_scale=18)


def draw_posterior_operator(axis):
    center = np.array([710.0, 353.0])
    rng = np.random.default_rng(4021)
    for angle in np.linspace(0, 2 * np.pi, 42, endpoint=False):
        inner = 118 + rng.uniform(-4, 6)
        outer = 160 + rng.uniform(-15, 16)
        start = center + inner * np.array([np.cos(angle), np.sin(angle)])
        end = center + outer * np.array([np.cos(angle), np.sin(angle)])
        axis.plot([start[0], end[0]], [start[1], end[1]], color=TEAL_MID, linewidth=1.1, alpha=0.16, zorder=1)
    axis.add_patch(Circle(center, 108, facecolor=blend(TEAL_LIGHT, amount=0.67), edgecolor="none", alpha=0.36, zorder=2))
    for radius, color, linewidth, alpha, linestyle in (
        (119, TEAL, 1.2, 1.0, (0, (3, 3))),
        (108, TEAL_LIGHT, 7.0, 0.55, "-"),
        (98, TEAL, 2.4, 0.95, "-"),
        (81, TEAL_MID, 2.0, 0.85, "-"),
        (66, TEAL, 1.3, 0.9, "-"),
    ):
        axis.add_patch(Circle(center, radius, facecolor="none", edgecolor=color, linewidth=linewidth, alpha=alpha, linestyle=linestyle, zorder=3))
    nodes = np.array(
        [[664, 400], [712, 433], [759, 400], [681, 360], [712, 376], [744, 360], [665, 322], [712, 322], [759, 322], [712, 285]]
    )
    links = ((0, 1), (1, 2), (0, 3), (2, 5), (3, 4), (4, 5), (3, 6), (4, 7), (5, 8), (6, 7), (7, 8), (7, 9))
    for first, second in links:
        axis.plot(nodes[[first, second], 0], nodes[[first, second], 1], color=INK, linewidth=2.5, zorder=5)
    axis.scatter(nodes[:, 0], nodes[:, 1], s=66, facecolor=blend(TEAL_LIGHT, amount=0.25), edgecolor=TEAL, linewidth=2.2, zorder=6)
    axis.text(710, 554, "ONE POSTERIOR", ha="center", color=TEAL, fontsize=17.5, weight="bold")
    axis.text(710, 523, "OPERATOR", ha="center", color=TEAL, fontsize=17.5, weight="bold")
    arrow(axis, (829, 353), (900, 353), AMBER, linewidth=3.0, mutation_scale=18)


def draw_shield(axis):
    for radius, alpha in ((116, 0.035), (98, 0.055), (82, 0.075)):
        axis.add_patch(Circle((1048, 505), radius, facecolor=AMBER, edgecolor="none", alpha=alpha, zorder=1))
    vertices = np.array([[1048, 618], [1127, 575], [1118, 451], [1048, 397], [978, 451], [969, 575]])
    shield = Polygon(vertices, closed=True, facecolor=blend(AMBER_LIGHT, amount=0.60), edgecolor=AMBER, linewidth=5.0, zorder=4)
    axis.add_patch(shield)
    add_shadow(shield, offset=(3, -4), alpha=0.13)
    axis.add_patch(Arc((1048, 521), 68, 82, theta1=0, theta2=180, edgecolor=AMBER, linewidth=8, zorder=6))
    rounded_box(axis, (1009, 450), 78, 69, AMBER, facecolor=AMBER, linewidth=1.0, radius=8, shadow=True, zorder=6)
    axis.add_patch(Circle((1048, 485), 9, facecolor=WHITE, edgecolor="none", zorder=7))
    axis.add_patch(Polygon([[1044, 480], [1052, 480], [1050, 461], [1046, 461]], facecolor=WHITE, edgecolor="none", zorder=7))


def draw_validation(axis):
    axis.text(1048, 681, "FROZEN WEIGHTS", ha="center", color=AMBER, fontsize=17, weight="bold")
    rounded_box(axis, (901, 42), 294, 622, AMBER, facecolor="none", linewidth=1.7, radius=18, linestyle=(0, (3, 3)), zorder=2)
    draw_shield(axis)
    for index, label in enumerate(("CALIBRATION", "BOUNDARY RECOVERY", "REPLICATION", "MCMC-DAGAR")):
        y = 320 - index * 78
        rounded_box(axis, (917, y), 268, 48, "#a67824", facecolor=WHITE, linewidth=1.15, radius=4, shadow=index == 0, zorder=5)
        rounded_box(axis, (930, y + 13), 20, 20, AMBER, facecolor=AMBER, linewidth=0.8, radius=3, zorder=6)
        axis.plot([935, 940, 947], [y + 23, y + 18, y + 29], color=WHITE, linewidth=1.7, solid_capstyle="round", zorder=7)
        axis.text(960, y + 24, label, ha="left", va="center", color=INK, fontsize=12.0, weight="bold", zorder=7)


def draw_gate(axis):
    x_left, x_right = 1232, 1350
    y_bottom, y_top = 295, 435
    axis.plot([x_left, x_left], [y_bottom, y_top], color=INK, linewidth=7, zorder=5)
    axis.plot([x_right, x_right], [y_bottom, y_top], color=INK, linewidth=7, zorder=5)
    for x in (x_left, x_right):
        axis.add_patch(Rectangle((x - 8, y_bottom - 7), 18, 8, facecolor=INK, edgecolor="none", zorder=5))
        axis.add_patch(Rectangle((x - 8, y_top - 1), 18, 8, facecolor=INK, edgecolor="none", zorder=5))
        axis.plot([x - 5, x + 6], [y_top + 9, y_top + 9], color=AMBER, linewidth=3, zorder=6)
    for x in np.linspace(x_left + 10, x_right - 10, 10):
        crown_y = y_top + 15 - 0.0063 * (x - (x_left + x_right) / 2) ** 2
        axis.plot([x, x], [y_bottom + 3, crown_y], color=INK, linewidth=2.0, zorder=4)
    curve_x = np.linspace(x_left + 4, x_right - 4, 120)
    curve_y = y_top + 15 - 0.0063 * (curve_x - (x_left + x_right) / 2) ** 2
    axis.plot(curve_x, curve_y, color=INK, linewidth=3.2, zorder=5)
    medallion = Circle((1291, 353), 23, facecolor=INK, edgecolor=AMBER, linewidth=2.5, zorder=8)
    axis.add_patch(medallion)
    add_shadow(medallion, offset=(2, -2), alpha=0.18)
    axis.plot([1280, 1289, 1304], [353, 342, 367], color=AMBER, linewidth=3.3, solid_capstyle="round", zorder=9)
    arrow(axis, (1194, 353), (1228, 353), CORAL, linewidth=3.0, mutation_scale=18)


def load_application_graph(spec):
    if not spec.gpkg_path.exists() or not spec.adjacency_path.exists():
        raise FileNotFoundError(f"Missing application input for {spec.label}")
    gdf = gpd.read_file(spec.gpkg_path)
    if len(gdf) != spec.n_areas:
        raise ValueError(f"{spec.label}: expected {spec.n_areas} areas, found {len(gdf)}")
    if gdf.crs is not None and gdf.crs.is_geographic:
        projected_crs = gdf.estimate_utm_crs()
        if projected_crs is not None:
            gdf = gdf.to_crs(projected_crs)
    centroids = gdf.geometry.centroid
    coordinates = np.column_stack([centroids.x.to_numpy(), centroids.y.to_numpy()])
    adjacency = pd.read_csv(spec.adjacency_path).to_numpy(dtype=float)
    if adjacency.shape != (spec.n_areas, spec.n_areas):
        raise ValueError(f"{spec.label}: adjacency shape {adjacency.shape} is inconsistent")
    covariate = pd.to_numeric(gdf[spec.covariate], errors="raise").to_numpy(dtype=float)
    return coordinates, adjacency, covariate


def scale_coordinates(coordinates, center, width, height):
    centered = coordinates - coordinates.mean(axis=0, keepdims=True)
    spans = np.ptp(centered, axis=0)
    spans[spans == 0] = 1.0
    scaled = centered * min(width / spans[0], height / spans[1])
    scaled[:, 0] += center[0]
    scaled[:, 1] += center[1]
    return scaled


def draw_application_graph(axis, spec):
    raw_coordinates, adjacency, covariate = load_application_graph(spec)
    coordinates = scale_coordinates(raw_coordinates, spec.center, *spec.box)
    edge_indices = np.argwhere(np.triu(adjacency, k=1) > 0)
    dissimilarities = np.abs(covariate[edge_indices[:, 0]] - covariate[edge_indices[:, 1]])
    highlighted = set(np.argsort(dissimilarities)[-max(5, len(edge_indices) // 28) :].tolist())
    for index, (first, second) in enumerate(edge_indices):
        is_boundary = index in highlighted
        axis.plot(
            coordinates[[first, second], 0],
            coordinates[[first, second], 1],
            color=CORAL if is_boundary else EDGE,
            linewidth=1.65 if is_boundary else 0.62,
            alpha=0.98 if is_boundary else 0.75,
            solid_capstyle="round",
            zorder=5 if is_boundary else 3,
        )
    axis.scatter(coordinates[:, 0], coordinates[:, 1], s=spec.node_size, facecolor=INK, edgecolor=WHITE, linewidth=0.25, zorder=6)
    label_x = 1680
    axis.text(label_x, spec.center[1] + 5, spec.label, color=INK, fontsize=14.5, weight="bold", va="center")
    sample_size_x = 1846
    axis.text(sample_size_x, spec.center[1] + 5, rf"$N = {spec.n_areas}$", color=INK, fontsize=13.5, va="center")


def draw_reuse(axis):
    axis.text(1652, 694, "NO RETRAINING", ha="center", color=CORAL, fontsize=17, weight="bold")
    axis.plot([1388, 1388], [272, 586], color=CORAL, linewidth=3.2, solid_capstyle="round", zorder=6)
    axis.plot([1354, 1388], [353, 353], color=CORAL, linewidth=3.2, zorder=6)
    for y in (585, 427, 267):
        arrow(axis, (1388, y), (1453, y), CORAL, linewidth=3.0, mutation_scale=18)
    for spec in APPLICATIONS:
        draw_application_graph(axis, spec)
    axis.plot([1498, 1498, 1922, 1922], [201, 184, 184, 201], color=CORAL, linewidth=3.0, solid_capstyle="round", zorder=6)
    arrow(axis, (1710, 184), (1710, 147), CORAL, linewidth=3.0, mutation_scale=18)
    rounded_box(axis, (1485, 42), 451, 105, CORAL, facecolor=WHITE, linewidth=2.2, radius=8, shadow=True, zorder=4)
    nodes = np.array([[1513, 94], [1536, 112], [1567, 93], [1533, 68]])
    for index, (first, second) in enumerate(((0, 1), (1, 2), (0, 3), (2, 3))):
        axis.plot(nodes[[first, second], 0], nodes[[first, second], 1], color=CORAL if index == 1 else INK, linewidth=2.5, zorder=6)
    axis.scatter(nodes[:, 0], nodes[:, 1], s=50, facecolor=INK, edgecolor=WHITE, linewidth=0.6, zorder=7)
    axis.text(1638, 106, "JOINT POSTERIOR", color=CORAL, fontsize=17.5, weight="bold", va="center", zorder=7)
    axis.text(1638, 72, "EDGE PROBABILITIES", color=INK, fontsize=16, va="center", zorder=7)


def draw_background(axis):
    axis.add_patch(
        Rectangle(
            (0, 0),
            CANVAS_WIDTH,
            CANVAS_HEIGHT,
            facecolor=WHITE,
            edgecolor="none",
            zorder=0,
        )
    )


def render(output_path):
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.color": INK, "axes.linewidth": 0})
    # The epsilon avoids a one-pixel truncation caused by binary float rounding.
    figure = plt.figure(
        figsize=((CANVAS_WIDTH + 0.01) / OUTPUT_DPI, (CANVAS_HEIGHT + 0.01) / OUTPUT_DPI),
        dpi=OUTPUT_DPI,
    )
    axis = figure.add_axes([0, 0, 1, 1])
    axis.set_xlim(0, CANVAS_WIDTH)
    axis.set_ylim(0, CANVAS_HEIGHT)
    axis.set_aspect("equal")
    axis.axis("off")
    draw_background(axis)
    draw_stage_headers(axis)
    draw_training_graphs(axis)
    draw_network_block(axis)
    draw_posterior_operator(axis)
    draw_validation(axis)
    draw_gate(axis)
    draw_reuse(axis)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=OUTPUT_DPI, facecolor=PAPER, edgecolor="none", metadata={"Software": "Matplotlib"})
    plt.close(figure)
    with Image.open(output_path) as rendered:
        rendered.convert("RGB").save(output_path, format="PNG", optimize=True)
    print(f"Saved: {output_path}")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH, help="Destination PNG path.")
    return parser.parse_args()


def main(output_path=None):
    render(OUTPUT_PATH if output_path is None else Path(output_path))


if __name__ == "__main__":
    args = parse_args()
    main(args.output)
