import cairosvg

def svg_to_png(input_svg, output_png, dpi=300, scale=1.0):
    """
    将 SVG 转换为 PNG
    
    参数:
        input_svg: 输入的 SVG 文件路径
        output_png: 输出的 PNG 文件路径
        dpi: 分辨率，默认 300
        scale: 缩放倍数，默认 1.0
    """
    cairosvg.svg2png(
        url=input_svg,
        write_to=output_png,
        dpi=dpi,
        scale=scale
    )
    print(f"转换完成：{output_png}（DPI={dpi}, scale={scale}）")


# 示例用法
if __name__ == "__main__":
    svg_to_png(
        input_svg="input.svg",
        output_png="output.png",
        dpi=300,    # 分辨率
        scale=2.0   # 放大2倍
    )