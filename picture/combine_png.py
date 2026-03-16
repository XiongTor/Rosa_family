#!/usr/bin/env python3
# Author: Tao Xiong
# Date: 2026-02-02
# Description: Merge all PDFs in a folder into a 4x6 grid image
# Python Version: 3.x

from pdf2image import convert_from_path
from PIL import Image
import os, glob

# ================= 配置区 (Configuration Section) =================

pdf_dir = "./"                  # PDF 所在目录 (PDF directory)
out_png = "Equidistants_combined.png"  # 输出文件名 (Output filename)

poppler_path = r"D:\python\poppler-25.12.0\Library\bin"

TARGET_WIDTH = 3000             # 每个子图宽度（像素） (Width of each sub-image in pixels)
COLS = 2
ROWS = 3
GAP = 50                        # 子图间距（像素） (Gap between sub-images in pixels)
BG_COLOR = (255, 255, 255)

# =========================================


def resize_to_width(img, target_w):
    w, h = img.size
    new_h = int(h * target_w / w)
    return img.resize((target_w, new_h), Image.LANCZOS)


# 1️⃣ 收集 PDF (Collect PDFs)
pdfs = sorted(glob.glob(os.path.join(pdf_dir, "*.pdf")))

if len(pdfs) == 0:
    raise RuntimeError("No PDF files found")

# 最多只取 24 个 (Take at most 24 files)
pdfs = pdfs[:COLS * ROWS]

# 2️⃣ PDF → Image → resize (Convert PDF to Image and resize)
images = []

for pdf in pdfs:
    imgs = convert_from_path(
        pdf,
        dpi=300,                 # 这里 DPI 不重要，后面统一 resize (DPI doesn't matter here, will resize uniformly later)
        poppler_path=poppler_path
    )
    img = resize_to_width(imgs[0], TARGET_WIDTH)
    images.append(img)

# 3️⃣ 计算单元格尺寸（统一高度） (Calculate cell dimensions - uniform height)
cell_w = max(img.size[0] for img in images)
cell_h = max(img.size[1] for img in images)

# 4️⃣ 计算画布尺寸 (Calculate canvas dimensions)
canvas_w = COLS * cell_w + (COLS - 1) * GAP
canvas_h = ROWS * cell_h + (ROWS - 1) * GAP

canvas = Image.new("RGB", (canvas_w, canvas_h), BG_COLOR)

# 5️⃣ 粘贴到网格 (Paste to grid)
for idx, img in enumerate(images):
    row = idx // COLS
    col = idx % COLS

    x = col * (cell_w + GAP) + (cell_w - img.size[0]) // 2
    y = row * (cell_h + GAP) + (cell_h - img.size[1]) // 2

    canvas.paste(img, (x, y))

# 6️⃣ 保存 (Save)
canvas.save(out_png, dpi=(300, 300))

print(f"Saved: {out_png}")
