#!/bin/bash
# Author: Tao Xiong
# Date: 2025-08-27
# Description: Record the script to build the phylogenetic tree by hybsuite
# ==== 主体代码开始 ====
hybsuite full_pipeline \
-input_list ./namelist.txt \
-eas_dir ./hybpiper_hybsuite \
-output_dir ./HybSuite_out_phylopypruner \
-t Reference_353.fasta \
-OI 1234567b \
-min_sample_coverage 0.1 \
-min_locus_coverage 0.1 \
-min_length 200 \
-skip_stage 012 \
-nt 15 \
-process 5