#!/usr/bin/env python3
"""Build bisulfite-converted FASTA + bowtie2 indices."""
import os, subprocess, sys

GENOME_FA = "data/hg19/hg19.fa"
OUT_DIR = "data/bowtie2_indices"

def bisulfite_convert_top(seq, methylated=False):
    """Convert top strand: C->T (unmethylated) or C->T except CpG (methylated)."""
    result = []
    seq = seq.upper()
    for i, base in enumerate(seq):
        if base == 'C':
            if methylated and i + 1 < len(seq) and seq[i + 1] == 'G':
                result.append('C')  # Keep CpG C
            else:
                result.append('T')  # Convert all other C
        else:
            result.append(base)
    return ''.join(result)

def write_fasta(seq, name, fh, width=80):
    fh.write(">" + name + "\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Get chromosome list from .fai
    with open(GENOME_FA + ".fai") as f:
        chroms = [(line.split("\t")[0], int(line.split("\t")[1])) for line in f]

    # Only main chromosomes (chr1-22, X, Y)
    main_chroms = [(c, l) for c, l in chroms
                   if c.startswith("chr") and "_" not in c
                   and (c.replace("chr","").isdigit() or c in ("chrX", "chrY"))]
    main_chroms.sort(key=lambda x: (x[0] not in ("chrX","chrY"),
                     int(x[0].replace("chr","")) if x[0].replace("chr","").isdigit() else 999))

    print("Processing " + str(len(main_chroms)) + " chromosomes")

    for suffix, methylated in [("converted_unmethylated", False),
                                ("converted_methylated", True)]:
        out_fasta = os.path.join(OUT_DIR, suffix + ".fa")
        if os.path.exists(out_fasta):
            print("  " + suffix + ".fa already exists, skipping")
            continue
        print("  Creating " + suffix + ".fa (methylated=" + str(methylated) + ")...")
        with open(out_fasta, 'w') as fh:
            for chrom, length in main_chroms:
                print("    " + chrom + " (" + str(length) + " bp)...")
                result = subprocess.run(
                    ["samtools", "faidx", GENOME_FA, chrom],
                    capture_output=True, text=True)
                lines = result.stdout.strip().split("\n")
                if len(lines) < 2:
                    continue
                seq = "".join(lines[1:]).upper()
                converted = bisulfite_convert_top(seq, methylated=methylated)
                write_fasta(converted, chrom, fh)
        print("  " + suffix + ".fa complete")

    # Build bowtie2 indices
    for name in ["unconverted", "converted_unmethylated", "converted_methylated"]:
        index_prefix = os.path.join(OUT_DIR, name)
        if os.path.exists(index_prefix + ".1.bt2"):
            print("  Index " + name + " exists, skipping")
            continue
        fasta_path = GENOME_FA if name == "unconverted" else os.path.join(OUT_DIR, name + ".fa")
        print("  Building bowtie2 index: " + name + "...")
        subprocess.run(["bowtie2-build", "--threads", "4", fasta_path, index_prefix], check=True)
        print("  Index " + name + " built")

    print("\nAll bowtie2 indices built!")

if __name__ == "__main__":
    main()
