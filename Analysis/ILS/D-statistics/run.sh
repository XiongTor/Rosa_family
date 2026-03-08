#!/bin/bash
# Author: Tao Xiong
# Date: 2026-01-16
# Description: Run all generated Dsuite scripts in parallel with nohup

echo "Starting parallel execution with nohup..."

# 遍历所有 Dsuite 文件夹并用 nohup 在后台运行
for folder in Dsuite_*/; do
    folder_name=$(basename "${folder}")
    script_file="${folder}${folder_name}.sh"

    if [ -f "${script_file}" ]; then
        echo "Launching: ${script_file}"
        nohup bash -c "cd ${folder} && bash ${folder_name}.sh" &
    fi
done

echo "All scripts launched in background with nohup!"
echo "Check logs in each folder: Dsuite_*/Dsuite_*.log"