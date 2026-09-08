"""Match DAGARBayes agreement-map geometry to the main CARBayes maps."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
CAR_RESULTS = ROOT / "Real Data Analysis" / "results_ABI_vs_CARBayes"
DAGAR_RESULTS = ROOT / "Real Data Analysis" / "results_ABI_vs_DAGARBayes"

APPLICATIONS = ("glasgow", "california", "south_korea")
LEGEND_HEIGHT_FRACTION = 0.75
HORIZONTAL_PADDING_PX = 80
MAP_LEGEND_GAP_PX = 120


def _content_bbox(image: Image.Image, tolerance: int = 8):
    pixels = np.asarray(image.convert("RGBA"))
    background = pixels[0, 0].astype(int)
    mask = np.max(np.abs(pixels.astype(int) - background), axis=2) > tolerance
    rows, columns = np.where(mask)
    if columns.size == 0:
        return (0, 0, image.width, image.height)
    return (
        int(columns.min()),
        int(rows.min()),
        int(columns.max()) + 1,
        int(rows.max()) + 1,
    )


def _column_clusters(mask, min_pixels: int = 10):
    columns = np.flatnonzero(mask.sum(axis=0) > min_pixels)
    if columns.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(columns) > 1)
    starts = [columns[0], *[columns[index + 1] for index in breaks]]
    ends = [*[columns[index] for index in breaks], columns[-1]]
    return [(int(start), int(end) + 1) for start, end in zip(starts, ends)]


def split_layout(image: Image.Image, tolerance: int = 8, merge_gap: int = 60):
    pixels = np.asarray(image.convert("RGBA"))
    background = pixels[0, 0].astype(int)
    mask = np.max(np.abs(pixels.astype(int) - background), axis=2) > tolerance
    clusters = _column_clusters(mask)
    if len(clusters) < 2:
        return image.crop(_content_bbox(image)), None

    legend_start, _ = clusters[-1]
    index = len(clusters) - 2
    while index >= 0 and legend_start - clusters[index][1] <= merge_gap:
        legend_start = clusters[index][0]
        index -= 1
    if index < 0:
        return image.crop(_content_bbox(image)), None

    map_image = image.crop((0, 0, legend_start, image.height))
    legend_image = image.crop((legend_start, 0, image.width, image.height))
    return (
        map_image.crop(_content_bbox(map_image)),
        legend_image.crop(_content_bbox(legend_image)),
    )


def _resize_to_height(image: Image.Image, target_height: int):
    if image.height == target_height:
        return image
    target_width = max(1, round(image.width * target_height / image.height))
    return image.resize((target_width, target_height), Image.Resampling.LANCZOS)


def standardize_map(reference_path: Path, target_path: Path):
    with Image.open(reference_path) as reference_raw, Image.open(target_path) as target_raw:
        reference = reference_raw.convert("RGBA")
        target = target_raw.convert("RGBA")

    reference_map, _ = split_layout(reference)
    target_map, target_legend = split_layout(target)
    if target_legend is None:
        raise RuntimeError(f"Could not isolate the legend in {target_path}")

    target_map = _resize_to_height(target_map, reference_map.height)
    target_legend = _resize_to_height(
        target_legend,
        round(LEGEND_HEIGHT_FRACTION * reference_map.height),
    )

    width = (
        2 * HORIZONTAL_PADDING_PX
        + target_map.width
        + MAP_LEGEND_GAP_PX
        + target_legend.width
    )
    height = max(target_map.height, target_legend.height)
    background = tuple(int(value) for value in np.asarray(reference)[0, 0])
    canvas = Image.new("RGBA", (width, height), background)

    map_x = HORIZONTAL_PADDING_PX
    map_y = (height - target_map.height) // 2
    legend_x = map_x + target_map.width + MAP_LEGEND_GAP_PX
    legend_y = (height - target_legend.height) // 2
    canvas.paste(target_map, (map_x, map_y), target_map)
    canvas.paste(target_legend, (legend_x, legend_y), target_legend)
    canvas.convert("RGB").save(target_path, dpi=(300, 300))

    print(
        f"Standardized {target_path.name}: canvas={canvas.size}, "
        f"map={target_map.size}, legend={target_legend.size}"
    )


def main():
    for application in APPLICATIONS:
        reference_path = CAR_RESULTS / application / f"{application}_boundary_agreement.png"
        target_path = DAGAR_RESULTS / application / f"{application}_boundary_agreement.png"
        if not reference_path.exists():
            raise FileNotFoundError(f"Missing reference map: {reference_path}")
        if not target_path.exists():
            raise FileNotFoundError(f"Missing DAGARBayes map: {target_path}")
        standardize_map(reference_path, target_path)


if __name__ == "__main__":
    main()
