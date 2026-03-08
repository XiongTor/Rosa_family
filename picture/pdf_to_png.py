#!/usr/bin/env python3
# Author: Tao Xiong
# Date: 2026-02-02
# Description:
# Python Version: 3.x
from pdf2image import convert_from_path
from PyPDF2 import PdfReader
import os, glob

TARGET_WIDTH = 12000
poppler_path = r"D:\python\poppler-25.12.0\Library\bin"

out_dir = "pngs"
os.makedirs(out_dir, exist_ok=True)   # 👈 关键这一行

for pdf in glob.glob("*.pdf"):
    base = os.path.splitext(pdf)[0]

    reader = PdfReader(pdf)
    page = reader.pages[0]
    width_inch = float(page.mediabox.width) / 72

    dpi = int(TARGET_WIDTH / width_inch)

    images = convert_from_path(
        pdf,
        dpi=dpi,
        poppler_path=poppler_path
    )

    images[0].save(os.path.join(out_dir, f"{base}.png"), "PNG")