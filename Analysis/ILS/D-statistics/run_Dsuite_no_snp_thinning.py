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
#     python run_Dsuite_biallelic.py [-f] alignment.fasta sets.txt msa2vcf_biallelic.py tree.txt Dsuite
# --
import sys
import itertools
import subprocess
import shutil
import os
import glob

def read_sample_groups(filename):
    samples_to_groups = {}
    groups_to_samples = {}
    extra_sample = None
    
    with open(filename, 'r') as f:
        for line in f:
            sample, group = line.strip().split()
            if group == 'xxx':
                continue  # Skip rows where the group is 'xxx'
            if group == 'Outgroup':
                extra_sample = sample
                continue  # Skip adding Outgroup to groups_to_samples
            samples_to_groups[sample] = group
            if group not in groups_to_samples:
                groups_to_samples[group] = set()
            groups_to_samples[group].add(sample)
    
    return samples_to_groups, groups_to_samples, extra_sample

def filter_fasta(input_fasta, output_fasta, samples_to_keep):
    with open(input_fasta, 'r') as in_file, open(output_fasta, 'w') as out_file:
        write_sequence = False
        for line in in_file:
            if line.startswith('>'):
                sample_name = line.strip().split()[0][1:]  # Remove '>' and get the first part of the header
                write_sequence = sample_name in samples_to_keep
                if write_sequence:
                    out_file.write(line)
            elif write_sequence:
                out_file.write(line)

def run_snp_sites(fasta_file, snp_sites_path_or_image):
    output_vcf = f"{fasta_file}_snps.vcf"
    current_dir = os.getcwd()  # Get the current working directory

    # Run snp-sites to generate the VCF
    if shutil.which(snp_sites_path_or_image):
        cmd = [snp_sites_path_or_image, '-v', '-o', output_vcf, fasta_file]
        subprocess.run(cmd, check=True)
    else:
        cmd = [
            'docker', 'run', '--platform', 'linux/amd64',
            '-v', f'{current_dir}:/data',
            snp_sites_path_or_image,
            'snp-sites', '-v', '-o', f"/data/{output_vcf}", f"/data/{fasta_file}"
        ]
        subprocess.run(cmd, check=True)

    # 🔧 修改：取消SNP稀疏化步骤
    # 原逻辑是每隔500bp只保留1个SNP，适用于全基因组数据；
    # 对于短片段多序列比对（如直系同源基因），这会导致信息损失，因此移除。

    # 保留原始输出，不进行 thinning
    return output_vcf
        

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

def create_filtered_sets_file(input_sets_file, samples_to_keep):
    base_name = os.path.splitext(input_sets_file)[0]
    filtered_sets_file = f"{base_name}_filtered.txt"
    
    with open(input_sets_file, 'r') as infile, open(filtered_sets_file, 'w') as outfile:
        for line in infile:
            sample, group = line.strip().split()
            if sample in samples_to_keep:
                outfile.write(f"{sample}\t{group}\n")
    
    return filtered_sets_file

def append_files(src_file, dest_file):
    with open(src_file, 'r') as infile, open(dest_file, 'a') as outfile:
        header = infile.readline()
        if os.stat(dest_file).st_size == 0:
            outfile.write(header)  # Write header only if dest_file is empty
        for line in infile:
            outfile.write(line)

def clean_up_files(patterns):
    for pattern in patterns:
        for file in glob.glob(pattern):
            os.remove(file)

def print_progress_bar(iteration, total, length=50):
    percent = ("{0:.1f}").format(100 * (iteration / float(total)))
    filled_length = int(length * iteration // total)
    bar = '█' * filled_length + '-' * (length - filled_length)
    sys.stdout.write(f'\rProgress: |{bar}| {percent}% Complete')
    sys.stdout.flush()

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python script.py [-f] input_fasta input_sets_file snp_sites_path_or_image tree_file dsuite_path")
        sys.exit(1)

    force_run = '-f' in sys.argv
    if force_run:
        sys.argv.remove('-f')

    input_fasta = sys.argv[1]
    input_sets_file = sys.argv[2]
    snp_sites_path_or_image = sys.argv[3]
    tree_file = sys.argv[4]
    dsuite_path = sys.argv[5]

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
    
    # Check if output files already exist
    if not force_run:
        existing_files = [tree_file_output, bbaa_file_output, dmin_file_output, log_file, error_file]
        for file in existing_files:
            if os.path.exists(file):
                print(f"Error: Output file '{file}' already exists. Use '-f' to force run and append to existing files.")
                sys.exit(1)

    for i, combo in enumerate(combinations_of_three):
        
        combo_str = ', '.join(combo)
        
        samples_to_keep = set()
        for group in combo:
            samples_to_keep.update(groups_to_samples[group])
        samples_to_keep.add(extra_sample)

        output_fasta = f"combination_{i+1}.fasta"
        filter_fasta(input_fasta, output_fasta, samples_to_keep)

        vcf_file = run_snp_sites(output_fasta, snp_sites_path_or_image)
        
        # Run Dsuite on the VCF file generated by snp-sites
        output_prefix = f"combination_{i+1}"
        filtered_sets_file = create_filtered_sets_file(input_sets_file, samples_to_keep)
        run_dsuite(vcf_file, output_prefix, tree_file, dsuite_path, filtered_sets_file, log_file, error_file, combo_str)
        
        # Append results to output files
        append_files(f"{output_prefix}_tree.txt", tree_file_output)
        append_files(f"{output_prefix}_BBAA.txt", bbaa_file_output)
        append_files(f"{output_prefix}_Dmin.txt", dmin_file_output)
        
        # Clean up intermediate files
        clean_up_files([output_fasta, vcf_file, f"{output_prefix}_tree.txt", f"{output_prefix}_BBAA.txt", f"{output_prefix}_Dmin.txt", filtered_sets_file])

        # Print progress bar
        print_progress_bar(i + 1, total_combinations)

    # Print final progress
    print_progress_bar(total_combinations, total_combinations)
    print("\nProcessing complete.")
