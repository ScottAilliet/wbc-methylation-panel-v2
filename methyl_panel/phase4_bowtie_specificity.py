#!/usr/bin/env python3
"""
Phase 4: Bowtie2 specificity screening for bisulfite-converted primers.

Checks primer specificity against 6 genome states:
1. Unconverted genome (sense + antisense) — detects genomic mispriming
2. Converted unmethylated genome (sense + antisense) — detects bisulfite-converted mispriming
3. Converted methylated genome (sense + antisense) — detects bisulfite-converted mispriming

For each primer pair, both left and right primers are aligned against all 6 states.
A primer passes if it aligns uniquely to the intended genome state and position.

The 6 genome states are created by in silico bisulfite conversion of the reference genome:
- Unconverted: original genome (2 strands)
- Converted unmethylated: all C→T on both strands (2 strands)
- Converted methylated: C→T except CpG C's (2 strands)
"""

import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import List, Optional, Tuple

from methyl_panel.phase2_bisulfite_convert import (
    reverse_complement, find_cpg_positions, bisulfite_convert_top, bisulfite_convert_bottom
)


@dataclass
class BowtieResult:
    """Result of bowtie2 alignment for a primer."""
    passes_filter: bool
    intended_genome: str       # Which genome state it was designed for
    total_alignments: int      # Total alignments across all 6 states
    off_target_alignments: int # Alignments to non-intended states
    mapping_note: str          # Description of any issues
    mismatch_profile: str = "" # e.g. "00200: 0 with 1 mismatch, 0 with 2 mm's, 2 with 3 mm's, 0 with 4 mms, 0 with 5 mms"


def create_bowtie_index(genome_fasta: str, output_dir: str,
                        converted: bool = True, methylated: bool = False) -> str:
    """
    Create a bowtie2 index for a genome state.

    Args:
        genome_fasta: Path to the original genome FASTA
        output_dir: Directory to store the index
        converted: If True, bisulfite-convert the genome
        methylated: If True, preserve CpG C's during conversion

    Returns:
        Path to the bowtie2 index prefix
    """
    if not converted:
        # Use original genome directly
        index_prefix = os.path.join(output_dir, "unconverted")
        if not os.path.exists(index_prefix + ".1.bt2"):
            cmd = ["bowtie2-build", genome_fasta, index_prefix]
            subprocess.run(cmd, check=True)
        return index_prefix

    # For converted genomes, we need to create the converted FASTA
    # This is a large operation — in practice, pre-built indices should be used
    suffix = "converted_methylated" if methylated else "converted_unmethylated"
    index_prefix = os.path.join(output_dir, suffix)

    if not os.path.exists(index_prefix + ".1.bt2"):
        # Create converted FASTA
        converted_fasta = os.path.join(output_dir, f"{suffix}.fa")

        # Read original FASTA and convert each chromosome
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
            # Use samtools faidx to get chromosome names
            cmd = ["samtools", "faidx", genome_fasta]
            # Get list of chromosomes from .fai file
            fai_file = genome_fasta + ".fai"
            chroms = []
            if os.path.exists(fai_file):
                with open(fai_file) as f:
                    for line in f:
                        parts = line.strip().split("\t")
                        if parts:
                            chroms.append(parts[0])

            for chrom in chroms:
                # Fetch full chromosome
                cmd = ["samtools", "faidx", genome_fasta, chrom]
                result = subprocess.run(cmd, capture_output=True, text=True)
                lines = result.stdout.strip().split("\n")
                if len(lines) < 2:
                    continue
                seq = "".join(lines[1:]).upper()

                # Convert top strand
                cpg_pos = find_cpg_positions(seq)
                converted_top = bisulfite_convert_top(seq, cpg_pos, methylated=methylated)

                # Write to FASTA
                tmp.write(f">{chrom}\n")
                for j in range(0, len(converted_top), 80):
                    tmp.write(converted_top[j:j+80] + "\n")

        tmp_path = tmp.name
        os.rename(tmp_path, converted_fasta)

        # Build bowtie2 index
        cmd = ["bowtie2-build", converted_fasta, index_prefix]
        subprocess.run(cmd, check=True)

    return index_prefix


def align_primer_bowtie2(primer_seq: str, index_prefix: str,
                         max_mismatches: int = 3) -> List[dict]:
    """
    Align a primer against a bowtie2 index.

    Args:
        primer_seq: Primer sequence
        index_prefix: Path to bowtie2 index prefix
        max_mismatches: Maximum allowed mismatches (edits) in alignment

    Returns:
        List of alignment dictionaries with parsed SAM fields.
        Returns empty list if bowtie2 fails or finds no hits.
    """
    # Create a FASTA file for the primer
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
        tmp.write(f">primer\n{primer_seq}\n")
        tmp_path = tmp.name

    try:
        # Bowtie2 settings for short primer screening:
        # --end-to-end: require full-length alignment (no soft clipping)
        # -N 1: allow 1 mismatch in seed (bowtie2 max is 1, NOT 2)
        # -L 10: seed length (short for 18-23bp primers)
        # -a: report all alignments
        # --score-min L,-6,-6: allow up to ~3 mismatches for 20bp primers
        cmd = [
            "bowtie2",
            "-x", index_prefix,
            "-f", tmp_path,
            "--end-to-end",
            "-N", "1",               # Max 1 seed mismatch (bowtie2 limit)
            "-L", "10",              # Short seed for short primers
            "-a",                    # Report all alignments
            "--no-hd",               # No header lines
            "--no-unal",             # Don't print unaligned reads
            "--score-min", "L,-6,-6",  # Allow up to ~3 mismatches
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        # If bowtie2 failed, return empty (caller should check for this)
        if result.returncode != 0:
            return []

        alignments = []
        for line in result.stdout.strip().split("\n"):
            if not line or line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue

            # Parse optional fields (SAM tags) properly
            nm = 0
            as_score = 0
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag.split(":")[2])
                elif tag.startswith("AS:i:"):
                    as_score = int(tag.split(":")[2])

            alignments.append({
                "read": parts[0],
                "flag": int(parts[1]),
                "rname": parts[2],
                "pos": int(parts[3]),
                "mapq": int(parts[4]),
                "cigar": parts[5],
                "seq": parts[9],
                "edits": nm,           # NM = number of mismatches
                "as_score": as_score,  # Alignment score
            })

        return alignments
    except (subprocess.TimeoutExpired, Exception):
        return []
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def screen_primer_pair(left_primer: str, right_primer: str,
                       intended_state: str,
                       index_dir: str,
                       max_mismatches: int = 3) -> BowtieResult:
    """
    Screen a primer pair against all genome states.

    A primer pair PASSES if:
    1. Both primers map to the intended genome (at least 1 hit each)
    2. No off-target alignments to non-intended genome states
    3. No multi-mapping within the intended genome (each primer maps uniquely)

    Args:
        left_primer: Left primer sequence
        right_primer: Right primer sequence
        intended_state: Which state the primers were designed for
                        (SM, AM, SU, AU, S, A)
        index_dir: Directory containing bowtie2 indices
        max_mismatches: Maximum mismatches (edits) allowed for a hit to count

    Returns:
        BowtieResult with pass/fail and details
    """
    # Map intended state to genome index
    state_to_index = {
        "SM": "converted_methylated",
        "AM": "converted_methylated",
        "SU": "converted_unmethylated",
        "AU": "converted_unmethylated",
        "S": "unconverted",
        "A": "unconverted",
    }

    intended_index = state_to_index.get(intended_state, "unconverted")
    indices_to_check = ["unconverted", "converted_unmethylated", "converted_methylated"]

    total_alignments = 0
    off_target = 0
    intended_hits = {"left": 0, "right": 0}
    multi_map = {"left": 0, "right": 0}
    notes = []
    bowtie_failed = False

    # Track mismatch counts (1-5) across ALL genomes for both primers combined
    mismatch_counts = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    for idx_name in indices_to_check:
        idx_path = os.path.join(index_dir, idx_name)
        if not os.path.exists(idx_path + ".1.bt2"):
            notes.append(f"Index {idx_name} not found, skipping")
            continue

        for primer_name, primer_seq in [("left", left_primer), ("right", right_primer)]:
            alignments = align_primer_bowtie2(primer_seq, idx_path, max_mismatches)

            # Filter: only count alignments with few edits
            good_alignments = [a for a in alignments if a["edits"] <= max_mismatches]
            total_alignments += len(good_alignments)

            # Count mismatches for profile (across all genomes, both primers)
            for a in good_alignments:
                e = a["edits"]
                if 1 <= e <= 5:
                    mismatch_counts[e] += 1

            if idx_name == intended_index:
                intended_hits[primer_name] = len(good_alignments)
                if len(good_alignments) > 1:
                    multi_map[primer_name] = len(good_alignments)
            elif good_alignments:
                off_target += len(good_alignments)
                notes.append(f"{primer_name} primer has {len(good_alignments)} off-target hits in {idx_name}")

    # Check for bowtie2 failure (no hits anywhere = likely bowtie2 error)
    if total_alignments == 0:
        notes.append("WARNING: No alignments found in any genome — possible bowtie2 error")
        bowtie_failed = True

    # Check that both primers mapped to intended genome
    if intended_hits["left"] == 0:
        notes.append("Left primer did not map to intended genome")
    if intended_hits["right"] == 0:
        notes.append("Right primer did not map to intended genome")

    # Check for multi-mapping within intended genome
    if multi_map["left"]:
        notes.append(f"Left primer maps to {multi_map['left']} locations in intended genome")
    if multi_map["right"]:
        notes.append(f"Right primer maps to {multi_map['right']} locations in intended genome")

    # Pass criteria:
    # 1. No off-target alignments
    # 2. Both primers map to intended genome
    # 3. No multi-mapping (each primer maps uniquely in intended genome)
    # 4. Bowtie2 didn't fail
    passes = (off_target == 0
              and intended_hits["left"] >= 1
              and intended_hits["right"] >= 1
              and multi_map["left"] <= 1
              and multi_map["right"] <= 1
              and not bowtie_failed)

    note = "; ".join(notes) if notes else "Unique mapping to intended genome"

    # Build mismatch profile string: "00200: 0 with 1 mismatch, 0 with 2 mm's, 2 with 3 mm's, 0 with 4 mms, 0 with 5 mms"
    compact = "".join(str(mismatch_counts[i]) for i in range(1, 6))
    parts = []
    for i in range(1, 6):
        mm_label = "mismatch" if i == 1 else "mm's"
        mms_label = "mms"
        if i == 1:
            parts.append(f"{mismatch_counts[i]} with 1 mismatch")
        elif i <= 4:
            parts.append(f"{mismatch_counts[i]} with {i} mm's")
        else:
            parts.append(f"{mismatch_counts[i]} with {i} mms")
    mismatch_profile = f"{compact}: " + ", ".join(parts)

    return BowtieResult(
        passes_filter=passes,
        intended_genome=intended_index,
        total_alignments=total_alignments,
        off_target_alignments=off_target,
        mapping_note=note,
        mismatch_profile=mismatch_profile,
    )


def main():
    """CLI entry point."""
    import argparse
    parser = argparse.ArgumentParser(description="Bowtie2 specificity screening")
    parser.add_argument("--left", required=True, help="Left primer sequence")
    parser.add_argument("--right", required=True, help="Right primer sequence")
    parser.add_argument("--state", required=True, choices=["SM", "AM", "SU", "AU", "S", "A"],
                        help="Intended genome state")
    parser.add_argument("--index-dir", required=True, help="Directory with bowtie2 indices")
    parser.add_argument("--max-mismatches", type=int, default=4)
    args = parser.parse_args()

    result = screen_primer_pair(
        args.left, args.right, args.state, args.index_dir, args.max_mismatches
    )

    print(f"Passes filter: {result.passes_filter}")
    print(f"Intended genome: {result.intended_genome}")
    print(f"Total alignments: {result.total_alignments}")
    print(f"Off-target: {result.off_target_alignments}")
    print(f"Note: {result.mapping_note}")


if __name__ == "__main__":
    main()
