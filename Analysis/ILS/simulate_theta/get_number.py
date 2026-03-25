import pdfplumber
import re
import os
import pandas as pd

def extract_data(folder_path):
    results = []
    
    # 1. 正则表达式准备
    # 提取文件名中的：样本名(Amy) 和 sim后的数字(0.001)
    # 格式假设：QT_样本名_..._sim数字.pdf
    file_pattern = re.compile(r'QT_(.*?)_.*_sim([\d\.]+)\.pdf')
    
    # 提取 PDF 内容中的 R2 数值
    # 匹配 R2 = 0.0462 或 R² = 0.0462，仅抓取数字部分
    r2_pattern = re.compile(r'R[²2]\s*=\s*([\d\.]+)')

    files = [f for f in os.listdir(folder_path) if f.endswith('.pdf')]
    
    for filename in files:
        # 解析文件名
        file_match = file_pattern.search(filename)
        if not file_match:
            continue
            
        sample_name = file_match.group(1)  # 第二部分: Amy
        sim_value = file_match.group(2)    # sim后的数字: 0.001
        
        r2_value = "Not Found"
        file_path = os.path.join(folder_path, filename)
        
        try:
            with pdfplumber.open(file_path) as pdf:
                for page in pdf.pages:
                    text = page.extract_text()
                    if text:
                        content_match = r2_pattern.search(text)
                        if content_match:
                            r2_value = content_match.group(1) # 只取数字部分
                            break # 找到第一个就跳出页面循环
            
            results.append({
                "Sim_Value": sim_value,
                "R2_Value": r2_value,
                "Sample_Name": sample_name
            })
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")

    # 2. 生成表格并排序（按 Sim_Value 数值排序更方便查看）
    df = pd.DataFrame(results)
    if not df.empty:
        # 转为浮点数排序后再转回字符串，保证 0.001 在 0.01 前面
        df['Sim_Value'] = pd.to_numeric(df['Sim_Value'])
        df = df.sort_values(by='Sim_Value').reset_index(drop=True)
        
        output_file = "extracted_phylogenetic_data.csv"
        df.to_csv(output_file, index=False, encoding='utf-8-sig')
        print(f"成功！提取了 {len(df)} 条数据，已保存至: {output_file}")
    else:
        print("未发现匹配的文件或数据。")

if __name__ == "__main__":
    extract_data('.')