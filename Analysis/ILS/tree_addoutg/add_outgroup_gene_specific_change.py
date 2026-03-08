#!/usr/bin/env python3
# Author: Tao Xiong
# Date: 2026-01-16
# Description: Add outgroup to gene trees with gene-specific branch lengths based on nearest ingroup relatives.
"""
add_outgroup_gene_specific_multi.py
Graft specified outgroup to gene trees missing any outgroup, estimating gene-specific branch lengths
based on nearest ingroup relatives. Supports multiple candidate outgroups, but only grafts one specified
target outgroup when none exist.

Dependencies:
    - biopython
    - matplotlib

Usage:
    python add_outgroup_gene_specific_multi.py \
        --input_trees trees_dir/ \
        --outgroups Zelkova_schneideriana Morus_indica Elaeagnus_angustifolia \
        --target_outgroup Zelkova_schneideriana \
        --output_dir grafted_trees/ \
        --sensitivity 0.2 \
        --histogram
"""

import os
import argparse
from Bio import Phylo
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(description="Add outgroup to gene trees with gene-specific branch lengths")

    ### 修改 1：允许输入多个外类群 + 一个指定插入的目标外类群
    parser.add_argument("--outgroups", nargs='+', required=True,
                        help="Names of possible outgroup taxa (space-separated)")
    parser.add_argument("--target_outgroup", required=True,
                        help="Name of the outgroup taxon to graft when none are present")

    parser.add_argument("--input_trees", required=True, help="Directory containing gene trees in Newick format")
    parser.add_argument("--output_dir", required=True, help="Directory to save grafted trees")
    parser.add_argument("--sensitivity", type=float, default=0.0,
                        help="Fractional range for sensitivity test (e.g., 0.2 for ±20%%)")
    parser.add_argument("--histogram", action="store_true", help="Generate histogram of branch lengths")
    return parser.parse_args()


def collect_tree_dict(tree_files):
    tree_dict = {}
    for f in tree_files:
        gene_name = os.path.basename(f)
        tree_dict[gene_name] = Phylo.read(f, "newick")
    return tree_dict


def get_branch_length_from_relatives(tree, outgroup_name, ingroup_taxa, reference_trees):
    lengths = []
    for ref_tree in reference_trees:
        ref_names = [t.name for t in ref_tree.get_terminals()]
        if all(t in ref_names for t in ingroup_taxa) and outgroup_name in ref_names:
            mrca = ref_tree.common_ancestor(ingroup_taxa)
            outgroup_terminal = ref_tree.find_any(name=outgroup_name)
            length = ref_tree.distance(outgroup_terminal, mrca)
            lengths.append(length)
    if lengths:
        return np.mean(lengths)
    else:
        fallback_lengths = []
        for ref_tree in reference_trees:
            names = [t.name for t in ref_tree.get_terminals()]
            if outgroup_name in names:
                outgroup_terminal = ref_tree.find_any(name=outgroup_name)
                path_length = ref_tree.distance(outgroup_terminal)
                fallback_lengths.append(path_length)
        return np.mean(fallback_lengths) if fallback_lengths else None


def graft_outgroup(tree, outgroup_name, branch_length):
    outgroup = Phylo.BaseTree.Clade(name=outgroup_name, branch_length=branch_length)
    new_root = Phylo.BaseTree.Clade()
    new_root.clades.append(deepcopy(tree.root))
    new_root.clades.append(outgroup)
    new_root.clades[0].branch_length = 0.0
    tree.root = new_root
    return tree


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    tree_files = [os.path.join(args.input_trees, f) for f in os.listdir(args.input_trees)
                  if f.endswith(".tre") or f.endswith(".nwk")]

    tree_dict = collect_tree_dict(tree_files)

    ### 修改 2：参考树中如果存在任意外类群则纳入
    trees_with_outgroup = [
        t for t in tree_dict.values()
        if any(og in [x.name for x in t.get_terminals()] for og in args.outgroups)
    ]

    grafted_lengths = []
    original_lengths = []
    summary = []

    for gene_name, tree in tree_dict.items():
        names = [t.name for t in tree.get_terminals()]
        present_outgroups = [og for og in args.outgroups if og in names]

        ### 修改 3：判断是否已有外类群
        if present_outgroups:
            used_og = present_outgroups[0]
            grafted_tree = tree
            applied_length = tree.distance(tree.find_any(name=used_og))
            original_lengths.append(applied_length)
            status = "existing"  # 已有外类群
        else:
            ingroup_taxa = names
            ### 修改 4：仅使用指定的 target_outgroup 进行 graft
            applied_length = get_branch_length_from_relatives(tree, args.target_outgroup, ingroup_taxa, trees_with_outgroup)
            if applied_length is None:
                print(f"⚠️ No reference found for {gene_name}, skipping grafting.")
                continue

            if args.sensitivity > 0:
                multiplier = 1 + np.random.uniform(-args.sensitivity, args.sensitivity)
                applied_length *= multiplier

            grafted_tree = graft_outgroup(tree, args.target_outgroup, applied_length)
            grafted_lengths.append(applied_length)
            status = "grafted"  # 新增外类群

        out_file = os.path.join(args.output_dir, gene_name)
        Phylo.write(grafted_tree, out_file, "newick")

        ### 修改 5：summary 文件中增加状态列
        summary.append(f"{gene_name}\t{applied_length:.6f}\t{status}")

    # Save summary
    summary_file = os.path.join(args.output_dir, "grafted_branch_lengths.txt")
    with open(summary_file, "w") as outsum:
        outsum.write("Gene\tApplied_Outgroup_Branch_Length\tStatus\n")
        outsum.write("\n".join(summary))

    # Plot histogram
    if args.histogram:
        plt.figure(figsize=(10, 6))
        bins = 20
        all_lengths = original_lengths + grafted_lengths
        if len(all_lengths) == 0:
            print("No data for histogram.")
        else:
            max_count = 1.05 * max(np.histogram(all_lengths, bins=bins)[0])

            if original_lengths:
                plt.hist(original_lengths, bins=bins, alpha=0.6, label='Original Outgroup',
                         color='green', edgecolor='black')

            if grafted_lengths:
                grafted_array = np.array(grafted_lengths)
                if args.sensitivity > 0:
                    lower = grafted_array / (1 + args.sensitivity)
                    upper = grafted_array / (1 - args.sensitivity)
                    for l, u in zip(lower, upper):
                        plt.fill_between([l, u], 0, max_count, color='skyblue', alpha=0.2)
                plt.hist(grafted_lengths, bins=bins, alpha=0.6, label='Grafted Outgroup',
                         color='skyblue', edgecolor='black')

            plt.xlabel("Outgroup Branch Length")
            plt.ylabel("Number of Genes")
            plt.title("Distribution of Outgroup Branch Lengths with Sensitivity")
            plt.legend()
            hist_file = os.path.join(args.output_dir, "outgroup_branch_lengths_histogram_sensitivity.png")
            plt.savefig(hist_file, dpi=300)
            plt.close()
            print(f"Sensitivity histogram saved to {hist_file}")

    print(f"Processing complete. Grafted trees saved in {args.output_dir}")
    print(f"Summary file saved to {summary_file}")


if __name__ == "__main__":
    main()
