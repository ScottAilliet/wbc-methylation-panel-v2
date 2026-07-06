#!/usr/bin/env python3
"""
Phase 6: Common SNP screening for primers.

Screens primers against dbSNP to identify common variants that could:
1. Disrupt primer binding (SNP in primer region)
2. Create false methylation calls (SNP creating/destroying a CpG)
3. Affect 3' end binding (SNP in last 5 nucleotides — most critical)

Uses a pre-downloaded dbSNP VCF file (common variants, MAF ≥ 1%).

Scoring: 2-digit score where:
- First digit: SNP in primer body (0=no, 1=yes)
- Second digit: SNP in 3' end last 5 nt (0=no, 1=yes)
Score 00 = clean, 10 = body SNP, 01 = 3' SNP, 11 = both
"""

import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class SNPResult:
    """Result of SNP screening."""
    common_variant_score: int   # 2-digit score (00, 10, 01, 11)
    has_snp_in_body: bool       # SNP anywhere in primer
    has_snp_in_3prime: bool     # SNP in last 5 nucleotides
    snp_count: int              # Total SNPs found
    snp_details: List[dict]     # Details of each SNP


def screen_primer_snp(primer_seq: str, chrom: str, genomic_start: int,
                      dbsnp_vcf: str, maf_threshold: float = 0.01,
                      primer_is_reverse: bool = False) -> SNPResult:
    """
    Screen a primer against dbSNP for common variants.

    Args:
        primer_seq: Primer sequence
        chrom: Chromosome name
        genomic_start: Genomic start position of the primer (1-based)
        dbsnp_vcf: Path to dbSNP VCF file (can be .gz)
        maf_threshold: Minimum MAF to consider a SNP "common"
        primer_is_reverse: If True, primer maps to reverse strand

    Returns:
        SNPResult with score and details
    """
    primer_len = len(primer_seq)
    genomic_end = genomic_start + primer_len - 1

    # Query dbSNP using tabix
    region = f"{chrom}:{genomic_start}-{genomic_end}"

    cmd = ["tabix", dbsnp_vcf, region]
    result = subprocess.run(cmd, capture_output=True, text=True)

    snps = []
    for line in result.stdout.strip().split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue

        snp_pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        info = parts[7]

        # Extract MAF from INFO field
        maf = 0.0
        for field in info.split(";"):
            if field.startswith("AF=") or field.startswith("MAF="):
                try:
                    maf = float(field.split("=")[1])
                except ValueError:
                    pass
            elif field.startswith("GMAF="):
                try:
                    maf = float(field.split("=")[1])
                except ValueError:
                    pass

        if maf >= maf_threshold:
            # Calculate position within primer
            pos_in_primer = snp_pos - genomic_start
            if primer_is_reverse:
                pos_in_primer = primer_len - 1 - pos_in_primer

            # Check if in 3' end (last 5 nucleotides)
            in_3prime = pos_in_primer >= primer_len - 5

            snps.append({
                "position": snp_pos,
                "pos_in_primer": pos_in_primer,
                "ref": ref,
                "alt": alt,
                "maf": maf,
                "in_3prime": in_3prime,
            })

    # Calculate score
    has_body = len(snps) > 0
    has_3prime = any(s["in_3prime"] for s in snps)

    # 2-digit score: first digit = body SNP, second digit = 3' SNP
    score = int(has_body) * 10 + int(has_3prime)

    return SNPResult(
        common_variant_score=score,
        has_snp_in_body=has_body,
        has_snp_in_3prime=has_3prime,
        snp_count=len(snps),
        snp_details=snps,
    )


def screen_primer_pair_snp(left_primer: str, right_primer: str,
                           chrom: str, left_genomic_start: int,
                           right_genomic_start: int,
                           dbsnp_vcf: str,
                           maf_threshold: float = 0.01) -> SNPResult:
    """
    Screen both primers of a pair against dbSNP.

    Returns the worst-case score across both primers.
    """
    left_result = screen_primer_snp(
        left_primer, chrom, left_genomic_start, dbsnp_vcf, maf_threshold,
        primer_is_reverse=False
    )
    right_result = screen_primer_snp(
        right_primer, chrom, right_genomic_start, dbsnp_vcf, maf_threshold,
        primer_is_reverse=True
    )

    # Combine: worst case
    combined_score = max(left_result.common_variant_score, right_result.common_variant_score)
    combined_snps = left_result.snp_details + right_result.snp_details

    return SNPResult(
        common_variant_score=combined_score,
        has_snp_in_body=left_result.has_snp_in_body or right_result.has_snp_in_body,
        has_snp_in_3prime=left_result.has_snp_in_3prime or right_result.has_snp_in_3prime,
        snp_count=len(combined_snps),
        snp_details=combined_snps,
    )


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Common SNP screening")
    parser.add_argument("--primer", required=True, help="Primer sequence")
    parser.add_argument("--chrom", required=True, help="Chromosome")
    parser.add_argument("--start", type=int, required=True, help="Genomic start (1-based)")
    parser.add_argument("--dbsnp", required=True, help="Path to dbSNP VCF")
    parser.add_argument("--maf", type=float, default=0.01, help="MAF threshold")
    args = parser.parse_args()

    result = screen_primer_snp(args.primer, args.chrom, args.start, args.dbsnp, args.maf)

    print(f"Score: {result.common_variant_score:02d}")
    print(f"SNP in body: {result.has_snp_in_body}")
    print(f"SNP in 3' end: {result.has_snp_in_3prime}")
    print(f"Total SNPs: {result.snp_count}")
    for s in result.snp_details:
        print(f"  pos={s['position']} ref={s['ref']} alt={s['alt']} MAF={s['maf']} 3'={s['in_3prime']}")


if __name__ == "__main__":
    main()
