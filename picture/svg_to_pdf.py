from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF

def svg_to_pdf(input_svg, output_pdf):
    """
    将 SVG 转换为 PDF
    
    参数:
        input_svg: 输入的 SVG 文件路径
        output_pdf: 输出的 PDF 文件路径
    """
    drawing = svg2rlg(input_svg)
    renderPDF.drawToFile(drawing, output_pdf)
    print(f"转换完成：{output_pdf}")


# 示例用法
if __name__ == "__main__":
    svg_to_pdf(
        input_svg="species_sets_no_geneflow_BBAA_f4ratio_new_order.svg",
        output_pdf="species_sets_no_geneflow_BBAA_f4ratio_new_order.pdf"
    )