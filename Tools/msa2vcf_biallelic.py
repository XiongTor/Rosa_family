#!/usr/bin/env python3
"""
msa2vcf_biallelic.py

从 MSA fasta 提取 biallelic SNP 并输出 VCF (uncompressed).
用法:
    python msa2vcf_biallelic.py alignment.fasta out.vcf --min_cov 0.8 --min_maf 0.01

注意:
 - 任何非 A/C/G/T (包括 N, -, ?) 被视为缺失。
 - 只输出 biallelic sites (exactly two alleles among non-missing).
 - REF = major allele, ALT = the other allele.
"""
import sys
import argparse
from collections import Counter
from Bio import SeqIO

VALID = set(['A','C','G','T'])

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("msa", help="Input MSA fasta (all sequences same length)")
    p.add_argument("vcf_out", help="Output VCF file")
    p.add_argument("--min_cov", type=float, default=1.0,
                   help="Minimum fraction of samples that must have a base (0-1). Default 1.0 (all samples)")
    p.add_argument("--min_maf", type=float, default=0.0,
                   help="Minimum minor allele frequency (as fraction of non-missing). Default 0.0")
    p.add_argument("--skip_masked", action='store_true',
                   help="If set, treat lowercase letters as masked (missing). Default: case-insensitive.")
    return p.parse_args()

def main():
    args = parse_args()
    records = list(SeqIO.parse(args.msa, "fasta"))
    if len(records) == 0:
        sys.exit("No sequences found in MSA")
    n = len(records)
    L = len(records[0].seq)
    # check equal length
    for r in records:
        if len(r.seq) != L:
            sys.exit("ERROR: sequences have different lengths")

    sample_names = [r.id for r in records]
    min_cov_count = int(round(args.min_cov * n))
    out = open(args.vcf_out, 'w')

    # write minimal VCF header
    out.write("##fileformat=VCFv4.2\n")
    out.write('##source=msa2vcf_biallelic\n')
    out.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"] + sample_names
    out.write("\t".join(header) + "\n")

    pos = 1
    for i in range(L):
        col = [str(r.seq[i]).upper() for r in records]
        # treat lowercase as same (we already made upper)
        bases = []
        idx2base = []
        for b in col:
            if b in VALID:
                bases.append(b)
                idx2base.append(b)
            else:
                bases.append(None)
                idx2base.append(None)
        non_missing = [b for b in bases if b is not None]
        if len(non_missing) < min_cov_count:
            pos += 1
            continue
        cnt = Counter(non_missing)
        alleles = list(cnt.keys())
        if len(alleles) != 2:
            pos += 1
            continue
        # minor allele frequency filter
        total_nonmiss = sum(cnt.values())
        freqs = sorted([(a, cnt[a]) for a in cnt], key=lambda x: -x[1])  # descending by count
        ref = freqs[0][0]
        alt = freqs[1][0]
        maf = freqs[1][1] / total_nonmiss
        if maf < args.min_maf:
            pos += 1
            continue
        # make genotype strings
        gts = []
        for b in idx2base:
            if b is None:
                gts.append("./.")
            elif b == ref:
                gts.append("0/0")
            elif b == alt:
                gts.append("1/1")
            else:
                # unexpected (shouldn't happen since len(alleles)==2), treat missing
                gts.append("./.")
        info = f"NS={total_nonmiss}"
        row = [ "chr1", str(pos), ".", ref, alt, ".", "PASS", info, "GT" ] + gts
        out.write("\t".join(row) + "\n")
        pos += 1

    out.close()
    print("Done. Wrote", args.vcf_out)

if __name__ == "__main__":
    main()
