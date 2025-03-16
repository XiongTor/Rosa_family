# !/usr/bin/env python
# Author: Tao Xiong
# Date: 2025.03.14
# encoding: utf-8

import requests
import pandas as pd

# 定义函数：逆地理编码查询
def reverse_geocoding(lat, lng, ak):
    url = "https://api.map.baidu.com/reverse_geocoding/v3/"
    params = {
        "ak": ak,
        "output": "json",
        "coordtype": "bd09ll",
        "extensions_poi": "1",
        "sort_strategy": "distance",
        "entire_poi": "1",
        "location": f"{lat},{lng}"
    }
    response = requests.get(url=url, params=params)
    if response.status_code == 200:
        data = response.json()
        if data["status"] == 0:
            result = data["result"]
            address = result["formatted_address"]
            return address
        else:
            return f"Error: {data['message']}"
    else:
        return "Error: 请求失败"

# 读取经纬度信息表格
def read_coordinates(file_path):
    try:
        df = pd.read_csv(file_path)  # 假设输入文件是CSV格式
        return df
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return None

# 批量查询并输出结果
def batch_query_and_output(file_path, ak, output_path):
    df = read_coordinates(file_path)
    if df is None:
        return
    
    # 确保表格中有'latitude'和'longitude'列
    if 'latitude' not in df.columns or 'longitude' not in df.columns:
        print("表格中必须包含'latitude'和'longitude'列")
        return

    # 创建新列存储查询结果
    df['address'] = df.apply(lambda row: reverse_geocoding(row['latitude'], row['longitude'], ak), axis=1)

    # 输出结果到新的CSV文件
    df.to_csv(output_path, index=False, encoding='utf-8-sig')
    print(f"结果已保存到文件: {output_path}")

# 主程序
if __name__ == "__main__":
    # 替换为你的百度地图AK
    ak = "amPjmvaQpoboUFtavVhwpdOlLjsMVd0e"
    
    # 输入文件路径（包含经纬度信息的CSV文件）
    input_file_path = "demo.csv"
    
    # 输出文件路径（保存查询结果的CSV文件）
    output_file_path = "result_addresses.csv"
    
    batch_query_and_output(input_file_path, ak, output_file_path)