#!/bin/bash
# Author: Tao Xiong
# Date: 2025-08-27
# Description: Record the script to build the phylogenetic tree by hybsuite
# ==== 主体代码开始 ====
hybsuite full_pipeline \
-input_list ./namelist.txt \
-eas_dir ./hybpiper_hybsuite \
-output_dir ./hybsuite_leng100_0.1 \
-t Reference_353.fasta \
-min_sample_coverage 0.1 \
-min_locus_coverage 0.1 \
-min_length 0 \
-skip_stage 01 \
-nt 12 \
-process 5





hybsuite full_pipeline \
-input_list ./namelist_tribe.txt \
-eas_dir ./hybpiper_hybsuite \
-output_dir ./HybSuite_out_tribe \
-t Reference_353.fasta \
-OI 1234567b \
-min_sample_coverage 0.1 \
-min_locus_coverage 0.1 \
-min_length 200 \
-skip_stage 01 \
-nt 12 \
-process 5


hybsuite full_pipeline \
-input_list ./namelist.txt \
-eas_dir ./hybpiper_assmble \
-output_dir ./hybsuite_out_orthofinder_2 \
-t orthofinder_v3_singel_rf.fasta \
-OI 1234567b \
-min_sample_coverage 0.1 \
-min_locus_coverage 0.1 \
-min_length 0 \
-skip_stage 01 \
-nt 12 \
-process 5