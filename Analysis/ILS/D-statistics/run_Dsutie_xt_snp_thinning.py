#!/usr/bin/env python3
# Author: Tao Xiong
# Date: 2025-10-14
# Description: 
# Python Version: 3.12.5
# -----------------------------------------
# 自动流程：
# 1️⃣ 从 sets.txt 读取样本群体；
# 2️⃣ 自动枚举三群组合（每次+Outgroup）；
# 3️⃣ 提取对应序列；
# 4️⃣ 调用 msa2vcf_biallelic.py 提取二态SNP并生成VCF；
# 5️⃣ 调用 Dsuite Dtrios；
# 6️⃣ 汇总结果；
# 7️⃣ 清理中间文件。

# 使用示例：
#     python run_Dsuite_biallelic.py [-f] alignment.fasta sets.txt msa2vcf_biallelic.py tree.txt Dsuite [--min_cov 0.8 --min_maf 0.01]
# -----------------------------------------
# ==== 主体代码开始 ====

import sys
import itertools
import subprocess
import shutil
import os
import glob
import argparse


# ================================================================
# 1. 读取 sets.txt（样本-群体映射）
# ================================================================
def read_sample_groups(filename):
    samples_to_groups = {}
    groups_to_samples = {}
    extra_sample = None
    
    with open(filename, 'r') as f:
        for line in f:
            sample, group = line.strip().split()
            if group == 'xxx':
                continue  # 跳过无效样本
            if group == 'Outgroup':
                extra_sample = sample
                continue  # 外群不计入普通群体
            samples_to_groups[sample] = group
            if group not in groups_to_samples:
                groups_to_samples[group] = set()
            groups_to_samples[group].add(sample)
    
    return samples_to_groups, groups_to_samples, extra_sample


# ================================================================
# 2. 过滤 FASTA，只保留需要的样本
# ================================================================
def filter_fasta(input_fasta, output_fasta, samples_to_keep):
    with open(input_fasta, 'r') as in_file, open(output_fasta, 'w') as out_file:
        write_sequence = False
        for line in in_file:
            if line.startswith('>'):
                sample_name = line.strip().split()[0][1:]  # 去掉 >
                write_sequence = sample_name in samples_to_keep
                if write_sequence:
                    out_file.write(line)
            elif write_sequence:
                out_file.write(line)


# ================================================================
# 3. 调用 msa2vcf_biallelic.py 提取二态性 SNP
# ================================================================
def run_msa2vcf_biallelic(fasta_file, msa2vcf_path, min_cov, min_maf):
    """
    调用 msa2vcf_biallelic.py 将 FASTA 转换为 VCF（仅保留二态性 SNP），
    然后对 SNP 进行稀疏化（每 500 bp 保留一个）。
    """
    output_vcf = f"{fasta_file}_biallelic.vcf"
    cmd = [
        "python", msa2vcf_path, fasta_file, output_vcf,
        "--min_cov", str(min_cov),
        "--min_maf", str(min_maf)
    ]
    subprocess.run(cmd, check=True)

    # SNP thinning：每500bp保留1个位点
    with open(output_vcf, 'r') as vcf_in:
        lines = vcf_in.readlines()

    with open(output_vcf, 'w') as vcf_out:
        previous_position = -500
        for line in lines:
            if line.startswith("#"):
                vcf_out.write(line)
            else:
                parts = line.split('\t')
                position = int(parts[1])
                if position >= previous_position + 500:
                    vcf_out.write(line)
                    previous_position = position

    return output_vcf


# ================================================================
# 4. 运行 Dsuite 分析
# ================================================================
def run_dsuite(vcf_file, output_prefix, tree_file, dsuite_path, input_names_file, log_file, error_file, combo):
    with open(log_file, 'a') as log, open(error_file, 'a') as err:
        log.write(f"\nRunning Dsuite for combination: {combo}\n")
        err.write(f"\nRunning Dsuite for combination: {combo}\n")
        log.flush()
        err.flush()
        cmd = [
            dsuite_path, 'Dtrios', '-k', '40', '-c',
            '-o', output_prefix, '-t', tree_file, vcf_file, input_names_file
        ]
        subprocess.run(cmd, stdout=log, stderr=err, check=True)


# ================================================================
# 5. 生成过滤后的 sets 文件
# ================================================================
def create_filtered_sets_file(input_sets_file, samples_to_keep):
    base_name = os.path.splitext(input_sets_file)[0]
    filtered_sets_file = f"{base_name}_filtered.txt"
    
    with open(input_sets_file, 'r') as infile, open(filtered_sets_file, 'w') as outfile:
        for line in infile:
            sample, group = line.strip().split()
            if sample in samples_to_keep:
                outfile.write(f"{sample}\t{group}\n")
    
    return filtered_sets_file


# ================================================================
# 6. 追加结果文件
# ================================================================
def append_files(src_file, dest_file):
    with open(src_file, 'r') as infile, open(dest_file, 'a') as outfile:
        header = infile.readline()
        if os.stat(dest_file).st_size == 0:
            outfile.write(header)  # 第一次写入表头
        for line in infile:
            outfile.write(line)


# ================================================================
# 7. 删除中间文件
# ================================================================
def clean_up_files(patterns):
    for pattern in patterns:
        for file in glob.glob(pattern):
            os.remove(file)


# ================================================================
# 8. 打印进度条
# ================================================================
def print_progress_bar(iteration, total, length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = '█' * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\rProgress: |{bar}| {percent}% Complete')
    sys.stdout.flush()


# ================================================================
# 9. 主程序入口
# ================================================================
if __name__ == "__main__":

    # 🔧修改点4：使用 argparse 支持 min_cov 和 min_maf 参数
    parser = argparse.ArgumentParser(
        description="自动化运行 Dsuite 分析（使用 msa2vcf_biallelic.py 生成二态SNP VCF）"
    )
    parser.add_argument("input_fasta", help="输入多序列比对文件（FASTA）")
    parser.add_argument("input_sets_file", help="样本-群体映射文件（两列：Sample Group）")
    parser.add_argument("msa2vcf_path", help="msa2vcf_biallelic.py 路径")
    parser.add_argument("tree_file", help="系统发育树文件（Newick）")
    parser.add_argument("dsuite_path", help="Dsuite 可执行路径")
    parser.add_argument("--min_cov", type=float, default=0.8, help="最小样本覆盖比例（默认0.8）")
    parser.add_argument("--min_maf", type=float, default=0.01, help="最小次等位频率（默认0.01）")
    parser.add_argument("-f", "--force", action="store_true", help="强制覆盖已有结果")

    args = parser.parse_args()

    input_fasta = args.input_fasta
    input_sets_file = args.input_sets_file
    msa2vcf_path = args.msa2vcf_path
    tree_file = args.tree_file
    dsuite_path = args.dsuite_path
    min_cov = args.min_cov
    min_maf = args.min_maf
    force_run = args.force

    samples_to_groups, groups_to_samples, extra_sample = read_sample_groups(input_sets_file)
    if extra_sample is None:
        print("Error: No Outgroup found in input_sets_file.")
        sys.exit(1)
    
    groups = list(groups_to_samples.keys())
    combinations_of_three = list(itertools.combinations(groups, 3))

    total_combinations = len(combinations_of_three)
    print(f"There are {total_combinations} total tests to run.")
    
    base_name = os.path.splitext(input_sets_file)[0]
    tree_file_output = f"{base_name}_tree.txt"
    bbaa_file_output = f"{base_name}_BBAA.txt"
    dmin_file_output = f"{base_name}_Dmin.txt"
    log_file = f"{base_name}_dsuite.log"
    error_file = f"{base_name}_dsuite.err"
    
    if not force_run:
        existing_files = [tree_file_output, bbaa_file_output, dmin_file_output, log_file, error_file]
        for file in existing_files:
            if os.path.exists(file):
                print(f"Error: Output file '{file}' already exists. Use -f to force run.")
                sys.exit(1)

    for i, combo in enumerate(combinations_of_three):
        combo_str = ', '.join(combo)
        samples_to_keep = set()
        for group in combo:
            samples_to_keep.update(groups_to_samples[group])
        samples_to_keep.add(extra_sample)

        output_fasta = f"combination_{i+1}.fasta"
        filter_fasta(input_fasta, output_fasta, samples_to_keep)

        vcf_file = run_msa2vcf_biallelic(output_fasta, msa2vcf_path, min_cov, min_maf)

        output_prefix = f"combination_{i+1}"
        filtered_sets_file = create_filtered_sets_file(input_sets_file, samples_to_keep)
        run_dsuite(vcf_file, output_prefix, tree_file, dsuite_path, filtered_sets_file, log_file, error_file, combo_str)
        
        append_files(f"{output_prefix}_tree.txt", tree_file_output)
        append_files(f"{output_prefix}_BBAA.txt", bbaa_file_output)
        append_files(f"{output_prefix}_Dmin.txt", dmin_file_output)
        
        clean_up_files([
            output_fasta, vcf_file,
            f"{output_prefix}_tree.txt",
            f"{output_prefix}_BBAA.txt",
            f"{output_prefix}_Dmin.txt",
            filtered_sets_file
        ])

        print_progress_bar(i + 1, total_combinations)

    print_progress_bar(total_combinations, total_combinations)
    print("\nProcessing complete.")
