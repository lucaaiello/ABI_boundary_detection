"""Normalize map legends and refresh the map assets used by the manuscript."""

from __future__ import annotations

from pathlib import Path
from shutil import copy2

import numpy as np
from PIL import Image

from standardize_dagar_agreement_maps import split_layout


ROOT = Path(__file__).resolve().parents[1]
REAL_DATA_DIR = ROOT / "Real Data Analysis"
EXPLORATORY_DIR = REAL_DATA_DIR / "results_exploratory_analysis"
CAR_DIR = REAL_DATA_DIR / "results_ABI_vs_CARBayes"
DAGAR_DIR = REAL_DATA_DIR / "results_ABI_vs_DAGARBayes"
MANUSCRIPT_IMAGE_DIR = ROOT / "Figure Generation" / "Images"

APPLICATIONS = ("glasgow", "california", "south_korea")
LEGEND_HEIGHT_FRACTION = 0.75
HORIZONTAL_PADDING_PX = 80
MAP_LEGEND_GAP_PX = 120

EXPLORATORY_FILES = (
    "glasgow_shr_map.png",
    "glasgow_income_deprivation_map.png",
    "california_lung_sir_map.png",
    "california_smoking_prevalence_map.png",
    "south_korea_lung_smr_map.png",
    "south_korea_smoking_prevalence_map.png",
)


def _resize_to_height(image: Image.Image, target_height: int):
    if image.height == target_height:
        return image
    target_width = max(1, round(image.width * target_height / image.height))
    return image.resize((target_width, target_height), Image.Resampling.LANCZOS)


def normalize_visible_legend(path: Path, reference_path: Path | None = None):
    with Image.open(path) as raw_image:
        image = raw_image.convert("RGBA")
    map_image, legend_image = split_layout(image)
    if legend_image is None:
        raise RuntimeError(f"Could not isolate the legend in {path}")

    if reference_path is None:
        canvas_height = image.height
    else:
        with Image.open(reference_path) as raw_reference:
            reference = raw_reference.convert("RGBA")
        reference_map, _ = split_layout(reference)
        map_image = _resize_to_height(map_image, reference_map.height)
        canvas_height = reference.height

    legend_height = round(LEGEND_HEIGHT_FRACTION * map_image.height)
    legend_image = _resize_to_height(legend_image, legend_height)

    canvas_width = (
        2 * HORIZONTAL_PADDING_PX
        + map_image.width
        + MAP_LEGEND_GAP_PX
        + legend_image.width
    )
    background = tuple(int(value) for value in np.asarray(image)[0, 0])
    canvas = Image.new("RGBA", (canvas_width, canvas_height), background)

    map_x = HORIZONTAL_PADDING_PX
    map_y = (canvas_height - map_image.height) // 2
    legend_x = map_x + map_image.width + MAP_LEGEND_GAP_PX
    legend_y = (canvas_height - legend_image.height) // 2
    canvas.paste(map_image, (map_x, map_y), map_image)
    canvas.paste(legend_image, (legend_x, legend_y), legend_image)
    canvas.convert("RGB").save(path, dpi=(300, 300))

    print(
        f"Normalized {path.name}: canvas={canvas.size}, map={map_image.size}, "
        f"legend={legend_image.size}, ratio={legend_image.height / map_image.height:.3f}"
    )


def copy_asset(source: Path, destination_name: str | None = None):
    destination = MANUSCRIPT_IMAGE_DIR / (destination_name or source.name)
    copy2(source, destination)
    print(f"Copied {source.relative_to(ROOT)} -> {destination.relative_to(ROOT)}")


def main():
    MANUSCRIPT_IMAGE_DIR.mkdir(parents=True, exist_ok=True)

    # Establish the main comparison maps first; all other panels reuse their
    # geographic map height while retaining their own aspect ratio.
    for application in APPLICATIONS:
        agreement = CAR_DIR / application / f"{application}_boundary_agreement.png"
        normalize_visible_legend(agreement)

    for application in APPLICATIONS:
        agreement = CAR_DIR / application / f"{application}_boundary_agreement.png"
        probability = CAR_DIR / application / f"{application}_boundary_probability_map.png"
        dagar_agreement = DAGAR_DIR / application / f"{application}_boundary_agreement.png"
        normalize_visible_legend(probability, reference_path=agreement)
        normalize_visible_legend(dagar_agreement, reference_path=agreement)

    exploratory_references = {
        "glasgow": CAR_DIR / "glasgow" / "glasgow_boundary_agreement.png",
        "california": CAR_DIR / "california" / "california_boundary_agreement.png",
        "south_korea": CAR_DIR / "south_korea" / "south_korea_boundary_agreement.png",
    }
    for filename in EXPLORATORY_FILES:
        application = next(name for name in APPLICATIONS if filename.startswith(name))
        normalize_visible_legend(
            EXPLORATORY_DIR / filename,
            reference_path=exploratory_references[application],
        )

    for filename in EXPLORATORY_FILES:
        copy_asset(EXPLORATORY_DIR / filename)

    for application in APPLICATIONS:
        copy_asset(CAR_DIR / application / f"{application}_boundary_probability_map.png")
        copy_asset(CAR_DIR / application / f"{application}_boundary_agreement.png")
        copy_asset(
            DAGAR_DIR / application / f"{application}_boundary_agreement.png",
            destination_name=f"{application}_boundary_agreement_dagarbayes.png",
        )


if __name__ == "__main__":
    main()
