#!/usr/bin/env python3
"""
Phase 3: Primer3-based primer design for bisulfite-converted DMR sequences.

Designs MSP (methylation-specific PCR) primers using Primer3 with
Primer3Plus-compatible settings. For each DMR block:

1. Fetch genomic sequence with flanking regions
2. Bisulfite convert to 6 strands
3. Design primers on the appropriate strand:
   - Methylated primers: from SM (sense methylated) or AM (antisense methylated)
   - Unmethylated primers: from SU (sense unmethylated) or AU (antisense unmethylated)
4. Apply bisulfite-specific constraints:
   - At least 2 CpGs per primer, 4 total per pair
   - Terminal CpG at 3' end for methylation discrimination
   - Tm 58-62°C (user-specified, required)
   - Product size 60-150 bp
"""

import os
import sys
import primer3
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

from methyl_panel.config import Primer3PlusConfig, PipelineConfig
from methyl_panel.phase2_bisulfite_convert import (
    convert_sequence, fetch_genomic_sequence, BisulfiteStrands, find_cpg_positions
)


@dataclass
class PrimerPair:
    """A designed primer pair with all properties."""
    # Identification
    assay_id: str
    seq_id: str
    cell_type_id: str
    template_used: str  # SM, AM, SU, AU

    # Primer sequences
    left_primer: str
    right_primer: str
    left_primer_display: str  # With CpG notation
    right_primer_display: str

    # Positions (0-based, relative to template)
    left_start: int
    left_len: int
    right_start: int
    right_len: int

    # Tm and properties
    left_tm: float
    right_tm: float
    product_size: int

    # CpG counts
    c_total_tail: int   # Total CpGs in both primer 3' tails
    c_total: int        # Total CpGs in amplicon
    left_c_total: int   # CpGs in left primer
    right_c_total: int  # CpGs in right primer
    left_c_tail: int    # CpGs in left primer 3' tail
    right_c_tail: int   # CpGs in right primer 3' tail

    # Mismatch scores (for methylation specificity)
    sense_meth_mismatch_score: int
    sense_unmeth_mismatch_score: int
    anti_meth_mismatch_score: int
    anti_unmeth_mismatch_score: int

    # Primer3 quality metrics
    left_gc_percent: float
    right_gc_percent: float
    left_self_any: float
    right_self_any: float
    left_self_end: float
    right_self_end: float
    pair_compl_any: float
    pair_compl_end: float
    left_end_stability: float
    right_end_stability: float
    penalty: float

    # Filter results (filled by later phases)
    bowtie_passes_filter: Optional[bool] = None
    bowtie_intended_genome: Optional[str] = None
    left_structure_mfe: Optional[float] = None
    right_structure_mfe: Optional[float] = None
    primer_dimer_prediction: Optional[str] = None
    primer_dimer_end_min_dg: Optional[float] = None
    common_variant_score: Optional[int] = None
    mapping_error_note: Optional[str] = None


def count_cpgs_by_position(cpg_positions: List[int], template_len: int,
                            template_name: str,
                            start: int, length: int) -> Tuple[int, int]:
    """
    Count CpG sites in a primer region using known genomic CpG positions.

    Bisulfite conversion is a substitution (C→T), not an indel, so CpG positions
    are identical between the genomic sequence and all bisulfite-converted
    templates. We use the known CpG positions from the genomic sequence rather
    than scanning the converted primer sequence, because:

    - In unmethylated templates (SU/AU), CpG C's become T's, so scanning for
      'CG' misses them. Scanning for 'TG' instead is wrong because non-CpG C's
      also become T's, creating many false 'TG' dinucleotides.
    - Position-based counting is template-independent and always correct.

    For sense strands (SM/SU), cpg_positions are top-strand 0-based positions
    and map directly to template coordinates.

    For antisense strands (AM/AU), the bottom strand is the reverse complement,
    so a top-strand CpG at position p maps to bottom-strand position len-p-2.

    Args:
        cpg_positions: 0-based positions of CpG C's on the top (sense) strand
        template_len: Length of the template sequence
        template_name: "SM", "AM", "SU", or "AU"
        start: 0-based start position of the primer in the template
        length: Length of the primer

    Returns:
        (total_cpg, tail_cpg) where tail_cpg = CpGs in last 5 nucleotides
    """
    # Convert top-strand CpG positions to this template's coordinate system
    if template_name in ("SM", "SU"):
        # Sense strand: positions map directly
        positions = cpg_positions
    else:
        # Antisense strand: reverse-complement mapping
        # A top-strand CpG at position p has its C at p and G at p+1.
        # On the bottom strand (reverse complement), the G becomes the
        # first base and the C becomes the second, at position len-p-2.
        positions = [template_len - p - 2 for p in cpg_positions]

    end = start + length  # exclusive

    # Count CpGs in the full primer region
    total = sum(1 for pos in positions if start <= pos < end)

    # Count CpGs in the 3' tail (last 5 nucleotides)
    tail_start = max(start, end - 5)
    tail_cpg = sum(1 for pos in positions if tail_start <= pos < end)

    return total, tail_cpg


def calculate_mismatch_score(primer_seq: str, template_strand: str,
                              primer_start: int) -> int:
    """
    Calculate mismatch score between primer and a template strand.
    Higher score = more mismatches = better discrimination.

    For MSP: a methylated primer should have many mismatches against
    the unmethylated template, and vice versa.
    """
    primer = primer_seq.upper()
    template = template_strand.upper()[primer_start:primer_start + len(primer)]

    if len(template) < len(primer):
        return len(primer)

    mismatches = sum(1 for p, t in zip(primer, template) if p != t)
    return mismatches


def design_primers_for_block(genomic_seq: str, config: Primer3PlusConfig,
                              pipeline_config: PipelineConfig,
                              assay_id: str, seq_id: str,
                              cell_type_id: str) -> List[PrimerPair]:
    """
    Design primers for a DMR block using Primer3.

    Tries both methylated (SM/AM) and unmethylated (SU/AU) templates.
    """
    # Bisulfite convert
    strands = convert_sequence(genomic_seq)

    # Validate Tm before running
    config.validate_tm()

    # Get Primer3 global args
    global_args = config.to_primer3_global_args()

    all_primers = []

    # Design on each of the 4 bisulfite-converted templates
    for template_name in ["SM", "AM", "SU", "AU"]:
        template_seq = strands.get_strand(template_name)

        # Skip if template is too short
        if len(template_seq) < 100:
            continue

        # Set up sequence args for Primer3
        seq_args = {
            "SEQUENCE_ID": f"{assay_id}_{template_name}",
            "SEQUENCE_TEMPLATE": template_seq,
        }

        # Run Primer3
        try:
            result = primer3.design_primers(seq_args, global_args)
        except Exception as e:
            continue

        n_returned = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n_returned == 0:
            continue

        # Process each returned primer pair
        for i in range(min(n_returned, pipeline_config.max_primers_per_block)):
            left_seq = result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
            right_seq = result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
            left_start, left_len = result.get(f"PRIMER_LEFT_{i}", (0, 0))
            right_start, right_len = result.get(f"PRIMER_RIGHT_{i}", (0, 0))

            if not left_seq or not right_seq:
                continue

            # Count CpGs using known genomic positions (template-independent)
            template_len = len(template_seq)
            left_cpg, left_cpg_tail = count_cpgs_by_position(
                strands.cpg_positions, template_len, template_name,
                left_start, left_len)
            right_cpg, right_cpg_tail = count_cpgs_by_position(
                strands.cpg_positions, template_len, template_name,
                right_start, right_len)

            # Skip if not enough CpGs (bisulfite-specific constraint)
            if (left_cpg + right_cpg) < pipeline_config.min_cpg_pair_total:
                continue
            if left_cpg < pipeline_config.min_cpg_per_primer and right_cpg < pipeline_config.min_cpg_per_primer:
                continue

            # Calculate mismatch scores against all 4 strands
            # For a methylated primer (from SM), it should mismatch SU
            sm = strands.sense_methylated
            su = strands.sense_unmethylated
            am = strands.antisense_methylated
            au = strands.antisense_unmethylated

            # Determine which strand the primer was designed from
            # and calculate mismatch against the opposite methylation state
            if template_name in ("SM", "AM"):
                # Methylated primer: check mismatch against unmethylated
                if template_name == "SM":
                    sm_mismatch = 0  # Perfect match to SM
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start) if template_name == "SM" else 0
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
                else:  # AM
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = 0
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
            else:  # SU or AU (unmethylated primers)
                if template_name == "SU":
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = 0
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start)
                    au_mismatch = calculate_mismatch_score(left_seq, au, left_start)
                else:  # AU
                    sm_mismatch = calculate_mismatch_score(left_seq, sm, left_start)
                    su_mismatch = calculate_mismatch_score(left_seq, su, left_start)
                    am_mismatch = calculate_mismatch_score(left_seq, am, left_start)
                    au_mismatch = 0

            # Create display sequences with CpG notation
            left_display = format_primer_display(left_seq)
            right_display = format_primer_display(right_seq)

            # Count total CpGs in amplicon using known genomic positions
            amp_start = left_start
            amp_end = right_start + right_len
            # Convert top-strand CpG positions to this template's coordinate system
            if template_name in ("SM", "SU"):
                amp_positions = strands.cpg_positions
            else:
                amp_positions = [template_len - p - 2 for p in strands.cpg_positions]
            amp_cpgs = sum(1 for pos in amp_positions if amp_start <= pos < amp_end)

            primer_pair = PrimerPair(
                assay_id=assay_id,
                seq_id=seq_id,
                cell_type_id=cell_type_id,
                template_used=template_name,
                left_primer=left_seq,
                right_primer=right_seq,
                left_primer_display=left_display,
                right_primer_display=right_display,
                left_start=left_start,
                left_len=left_len,
                right_start=right_start,
                right_len=right_len,
                left_tm=result.get(f"PRIMER_LEFT_{i}_TM", 0),
                right_tm=result.get(f"PRIMER_RIGHT_{i}_TM", 0),
                product_size=result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0),
                c_total_tail=left_cpg_tail + right_cpg_tail,
                c_total=amp_cpgs,
                left_c_total=left_cpg,
                right_c_total=right_cpg,
                left_c_tail=left_cpg_tail,
                right_c_tail=right_cpg_tail,
                sense_meth_mismatch_score=sm_mismatch,
                sense_unmeth_mismatch_score=su_mismatch,
                anti_meth_mismatch_score=am_mismatch,
                anti_unmeth_mismatch_score=au_mismatch,
                left_gc_percent=result.get(f"PRIMER_LEFT_{i}_GC_PERCENT", 0),
                right_gc_percent=result.get(f"PRIMER_RIGHT_{i}_GC_PERCENT", 0),
                left_self_any=result.get(f"PRIMER_LEFT_{i}_SELF_ANY", 0),
                right_self_any=result.get(f"PRIMER_RIGHT_{i}_SELF_ANY", 0),
                left_self_end=result.get(f"PRIMER_LEFT_{i}_SELF_END", 0),
                right_self_end=result.get(f"PRIMER_RIGHT_{i}_SELF_END", 0),
                pair_compl_any=result.get(f"PRIMER_PAIR_{i}_COMPL_ANY", 0),
                pair_compl_end=result.get(f"PRIMER_PAIR_{i}_COMPL_END", 0),
                left_end_stability=result.get(f"PRIMER_LEFT_{i}_END_STABILITY", 0),
                right_end_stability=result.get(f"PRIMER_RIGHT_{i}_END_STABILITY", 0),
                penalty=result.get(f"PRIMER_PAIR_{i}_PENALTY", 0),
            )

            all_primers.append(primer_pair)

    return all_primers


def format_primer_display(seq: str) -> str:
    """Format primer sequence with CpG notation (lowercase y for C in CpG)."""
    result = []
    seq_upper = seq.upper()
    for i, base in enumerate(seq_upper):
        if base == 'C' and i + 1 < len(seq_upper) and seq_upper[i + 1] == 'G':
            result.append('y')  # CpG C
        elif base == 'G' and i > 0 and seq_upper[i - 1] == 'C':
            result.append('g')  # CpG G (already shown by preceding y)
        else:
            result.append(base)
    return ''.join(result)


def main():
    """CLI entry point for testing."""
    import argparse
    parser = argparse.ArgumentParser(description="Design primers for DMR blocks")
    parser.add_argument("--genome", default="/workspace/wgbs_tools/references/hg19/hg19.fa.gz",
                        help="Path to genome FASTA")
    parser.add_argument("--region", help="Genomic region (chr:start-end)")
    parser.add_argument("--settings", help="Primer3Plus settings file")
    parser.add_argument("--min-tm", type=float, default=58.0)
    parser.add_argument("--opt-tm", type=float, default=60.0)
    parser.add_argument("--max-tm", type=float, default=62.0)
    args = parser.parse_args()

    # Load config
    if args.settings:
        config = Primer3PlusConfig.from_primer3plus_file(args.settings)
    else:
        config = Primer3PlusConfig()

    # Set Tm (required)
    config.primer_min_tm = args.min_tm
    config.primer_opt_tm = args.opt_tm
    config.primer_max_tm = args.max_tm

    pipeline_config = PipelineConfig()

    # Fetch sequence
    if args.region:
        chrom, coords = args.region.split(":")
        start, end = coords.split("-")
        seq = fetch_genomic_sequence(args.genome, chrom, int(start), int(end), flank=50)
    else:
        print("Please provide --region")
        sys.exit(1)

    print(f"Sequence length: {len(seq)}")
    print(f"Tm range: {config.primer_min_tm}-{config.primer_max_tm}°C")

    # Design primers
    primers = design_primers_for_block(
        seq, config, pipeline_config, "test_001", "TEST_0001", "TEST"
    )

    print(f"\nDesigned {len(primers)} primer pairs:")
    for p in primers[:5]:
        print(f"\n  {p.assay_id} (template: {p.template_used})")
        print(f"    Left:  {p.left_primer} (Tm={p.left_tm:.1f}°C, {p.left_c_total}CpGs)")
        print(f"    Right: {p.right_primer} (Tm={p.right_tm:.1f}°C, {p.right_c_total}CpGs)")
        print(f"    Product: {p.product_size}bp, Total CpGs: {p.c_total}")
        print(f"    Penalty: {p.penalty:.2f}")


if __name__ == "__main__":
    main()
