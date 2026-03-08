#!/usr/bin/env python3
# Author: Tao Xiong
# Date: 2026-03-04
# Description: Convert PDF files to PNG images with target width
# Usage: python pdf_to_png.py [--input DIR] [--output DIR] [--width INT] [--pages all|first]
# Dependencies: pdf2image, PyPDF2, poppler
# Python Version: 3.x

import os
import glob
import argparse
from pdf2image import convert_from_path
from PyPDF2 import PdfReader

POPPLER_PATH = r"D:\python\poppler-25.12.0\Library\bin"

def pdf_to_png(input_dir=".", output_dir="pngs", target_width=12000, pages="all"):
    os.makedirs(output_dir, exist_ok=True)

    pdf_files = glob.glob(os.path.join(input_dir, "*.pdf"))
    if not pdf_files:
        print(f"No PDF files found in: {input_dir}")
        return

    for pdf in pdf_files:
        base = os.path.splitext(os.path.basename(pdf))[0]
        print(f"Processing: {pdf}")

        reader = PdfReader(pdf)
        page = reader.pages[0]
        width_inch = float(page.mediabox.width) / 72
        dpi = int(target_width / width_inch)

        images = convert_from_path(pdf, dpi=dpi, poppler_path=POPPLER_PATH)

        if pages == "first":
            images = images[:1]

        for i, img in enumerate(images):
            suffix = f"_p{i+1}" if len(images) > 1 else ""
            out_path = os.path.join(output_dir, f"{base}{suffix}.png")
            img.save(out_path, "PNG")
            print(f"  Saved: {out_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert PDF files to PNG images")
    parser.add_argument("--input",  default=".",     help="Input directory (default: current dir)")
    parser.add_argument("--output", default="pngs",  help="Output directory (default: pngs/)")
    parser.add_argument("--width",  type=int, default=12000, help="Target width in pixels (default: 12000)")
    parser.add_argument("--pages",  choices=["all", "first"], default="all", help="Which pages to convert (default: all)")
    args = parser.parse_args()

    pdf_to_png(args.input, args.output, args.width, args.pages)
