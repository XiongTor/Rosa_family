#!/usr/bin/env python3
"""
Concatenate images into a fixed grid layout.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

from PIL import Image


RGBColor = tuple[int, int, int]


def parse_layout(layout: str) -> tuple[int, int]:
    try:
        rows_text, cols_text = layout.lower().split("x", maxsplit=1)
        rows = int(rows_text)
        cols = int(cols_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid layout '{layout}'. Use the form ROWSxCOLS, for example 2x3."
        ) from exc

    if rows <= 0 or cols <= 0:
        raise argparse.ArgumentTypeError("Layout rows and columns must be positive integers.")

    return rows, cols


def positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("Value must be a positive integer.")
    return number


def non_negative_int(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Value must be a non-negative integer.")
    return number


def parse_background(value: str) -> RGBColor:
    normalized = value.strip().lstrip("#")

    if len(normalized) != 6:
        raise argparse.ArgumentTypeError(
            "Background color must be a 6-digit hex value such as FFFFFF or #FFFFFF."
        )

    try:
        rgb = tuple(int(normalized[index : index + 2], 16) for index in range(0, 6, 2))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Background color must be a valid hexadecimal RGB value."
        ) from exc

    return rgb


def iter_image_paths(paths: Iterable[str]) -> list[Path]:
    image_paths = [Path(path) for path in paths]

    missing = [str(path) for path in image_paths if not path.is_file()]
    if missing:
        missing_text = "\n".join(f"  - {path}" for path in missing)
        raise FileNotFoundError(f"The following input files do not exist:\n{missing_text}")

    return image_paths


def inspect_grid(image_paths: list[Path], rows: int, cols: int) -> tuple[int, int]:
    expected = rows * cols
    if len(image_paths) != expected:
        raise ValueError(f"Layout {rows}x{cols} requires {expected} images, got {len(image_paths)}.")

    max_width = 0
    max_height = 0

    for path in image_paths:
        with Image.open(path) as img:
            max_width = max(max_width, img.width)
            max_height = max(max_height, img.height)

    return max_width, max_height


def normalize_image(image: Image.Image, background: RGBColor) -> Image.Image:
    if image.mode in {"RGBA", "LA"} or (image.mode == "P" and "transparency" in image.info):
        rgba_image = image.convert("RGBA")
        canvas = Image.new("RGBA", rgba_image.size, background + (255,))
        return Image.alpha_composite(canvas, rgba_image).convert("RGB")

    return image.convert("RGB")


def concat_images(
    image_paths: list[Path],
    output_path: Path,
    rows: int,
    cols: int,
    dpi: int = 300,
    gap: int = 0,
    background: RGBColor = (255, 255, 255),
) -> None:
    cell_width, cell_height = inspect_grid(image_paths, rows, cols)

    canvas_width = cell_width * cols + gap * (cols - 1)
    canvas_height = cell_height * rows + gap * (rows - 1)
    canvas = Image.new("RGB", (canvas_width, canvas_height), background)

    for index, path in enumerate(image_paths):
        row = index // cols
        col = index % cols

        with Image.open(path) as img:
            normalized = normalize_image(img, background)
            x = col * (cell_width + gap) + (cell_width - normalized.width) // 2
            y = row * (cell_height + gap) + (cell_height - normalized.height) // 2
            canvas.paste(normalized, (x, y))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output_path, dpi=(dpi, dpi))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Concatenate multiple images into a grid.")
    parser.add_argument("images", nargs="+", help="Input image paths in placement order.")
    parser.add_argument("-o", "--output", required=True, help="Output image path.")
    parser.add_argument(
        "-l",
        "--layout",
        required=True,
        help="Grid layout in the form ROWSxCOLS, for example 2x2 or 1x4.",
    )
    parser.add_argument(
        "-d",
        "--dpi",
        type=positive_int,
        default=300,
        help="Output DPI. Default: 300.",
    )
    parser.add_argument(
        "-g",
        "--gap",
        type=non_negative_int,
        default=0,
        help="Gap between grid cells in pixels. Default: 0.",
    )
    parser.add_argument(
        "-b",
        "--background",
        type=parse_background,
        default=(255, 255, 255),
        help="Background color as hex RGB, for example FFFFFF. Default: FFFFFF.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        rows, cols = parse_layout(args.layout)
        image_paths = iter_image_paths(args.images)
        output_path = Path(args.output)

        concat_images(
            image_paths=image_paths,
            output_path=output_path,
            rows=rows,
            cols=cols,
            dpi=args.dpi,
            gap=args.gap,
            background=args.background,
        )
    except (argparse.ArgumentTypeError, FileNotFoundError, OSError, ValueError) as exc:
        parser.exit(status=1, message=f"Error: {exc}\n")

    print(f"Saved: {output_path} (layout: {rows}x{cols}, dpi: {args.dpi})")


if __name__ == "__main__":
    main()
