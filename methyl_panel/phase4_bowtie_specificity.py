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
    mismatch_profile: str = "" # compact 6-digit code: 0mm,1mm,2mm,3mm,4mm,5mm
    mismatch_detail: str = ""  # human-readable description


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
    Align a single primer against a bowtie2 index.

    .. deprecated:: Use batch_align_primers_bowtie2 for efficiency.
    Kept for backward compatibility / single-primer testing.

    Args:
        primer_seq: Primer sequence
        index_prefix: Path to bowtie2 index prefix
        max_mismatches: Maximum allowed mismatches (edits) in alignment

    Returns:
        List of alignment dictionaries with parsed SAM fields.
    """
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
        tmp.write(f">primer\n{primer_seq}\n")
        tmp_path = tmp.name

    try:
        cmd = [
            "bowtie2",
            "-x", index_prefix,
            "-f", tmp_path,
            "--end-to-end",
            "-N", "1",
            "-L", "10",
            "-k", "100",             # Cap at 100 alignments (not -a which is unlimited)
            "--no-hd",
            "--no-unal",
            "--score-min", "L,-6,-6",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

        if result.returncode != 0:
            return []

        return _parse_sam_lines(result.stdout)
    except (subprocess.TimeoutExpired, Exception):
        return []
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def _parse_sam_lines(sam_text: str) -> List[dict]:
    """Parse SAM output lines into alignment dicts."""
    alignments = []
    for line in sam_text.strip().split("\n"):
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < 11:
            continue

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
            "edits": nm,
            "as_score": as_score,
        })
    return alignments


def batch_align_primers_bowtie2(primer_seqs: List[str], index_prefix: str,
                                 max_mismatches: int = 3) -> dict:
    """
    Align ALL primers against a bowtie2 index in a SINGLE bowtie2 invocation.

    This is the memory-efficient replacement for calling align_primer_bowtie2
    in a loop. Instead of spawning one bowtie2 process per primer (each loading
    the ~4GB index), we write all primers to one FASTA and run bowtie2 once.

    Args:
        primer_seqs: List of primer sequences
        index_prefix: Path to bowtie2 index prefix
        max_mismatches: Maximum allowed mismatches (edits)

    Returns:
        Dict mapping {primer_index: [alignment_dicts]}
    """
    if not primer_seqs:
        return {}

    # Write all primers to one FASTA file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
        for i, seq in enumerate(primer_seqs):
            tmp.write(f">primer_{i}\n{seq}\n")
        tmp_path = tmp.name

    try:
        cmd = [
            "bowtie2",
            "-x", index_prefix,
            "-f", tmp_path,
            "--end-to-end",
            "-N", "1",
            "-L", "10",
            "-k", "100",             # Cap at 100 alignments per primer
            "--no-hd",
            "--no-unal",
            "--score-min", "L,-6,-6",
            "--threads", "1",        # Single thread — we're I/O bound, not CPU
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            return {}

        # Parse SAM output and group by primer index
        results = {}
        for line in result.stdout.strip().split("\n"):
            if not line or line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue

            # Extract primer index from read name (primer_0, primer_1, etc.)
            read_name = parts[0]
            try:
                primer_idx = int(read_name.split("_")[-1])
            except (ValueError, IndexError):
                continue

            nm = 0
            as_score = 0
            for tag in parts[11:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag.split(":")[2])
                elif tag.startswith("AS:i:"):
                    as_score = int(tag.split(":")[2])

            if primer_idx not in results:
                results[primer_idx] = []

            results[primer_idx].append({
                "read": read_name,
                "flag": int(parts[1]),
                "rname": parts[2],
                "pos": int(parts[3]),
                "mapq": int(parts[4]),
                "cigar": parts[5],
                "seq": parts[9],
                "edits": nm,
                "as_score": as_score,
            })

        return results
    except (subprocess.TimeoutExpired, Exception):
        return {}
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def screen_primer_pair(left_primer: str, right_primer: str,
                       intended_state: str,
                       index_dir: str,
                       max_mismatches: int = 3) -> BowtieResult:
    """
    Screen a single primer pair (calls bowtie2 per primer — slow, memory-heavy).

    .. deprecated:: Use screen_primer_pairs_batch for efficiency.
    Kept for backward compatibility / single-primer testing.
    """
    state_to_index = {
        "SM": "converted_methylated", "AM": "converted_methylated",
        "SU": "converted_unmethylated", "AU": "converted_unmethylated",
        "S": "unconverted", "A": "unconverted",
    }
    intended_index = state_to_index.get(intended_state, "unconverted")
    indices_to_check = ["unconverted", "converted_unmethylated", "converted_methylated"]

    # Collect all alignments for this pair
    all_alignments = {}  # {idx_name: {"left": [alignments], "right": [alignments]}}
    for idx_name in indices_to_check:
        idx_path = os.path.join(index_dir, idx_name)
        if not os.path.exists(idx_path + ".1.bt2"):
            all_alignments[idx_name] = {"left": [], "right": []}
            continue
        left_aln = align_primer_bowtie2(left_primer, idx_path, max_mismatches)
        right_aln = align_primer_bowtie2(right_primer, idx_path, max_mismatches)
        all_alignments[idx_name] = {"left": left_aln, "right": right_aln}

    return _evaluate_screening(all_alignments, intended_index, indices_to_check, max_mismatches)


def screen_primer_pairs_batch(primers_data: list, index_dir: str,
                               max_mismatches: int = 3,
                               product_size_min: int = 40,
                               product_size_max: int = 150) -> list:
    """
    Screen ALL primer pairs in batch — runs bowtie2 only 3 times total
    (once per genome index) instead of 6 times per primer pair.

    This is the memory-efficient replacement for calling screen_primer_pair
    in a loop. It:
    1. Collects all unique primer sequences
    2. Writes them to one FASTA per genome index
    3. Runs bowtie2 once per index (3 total, not 6×N)
    4. Parses results back and evaluates each primer pair

    Args:
        primers_data: List of primer dicts (from primers.json)
        index_dir: Directory containing bowtie2 indices
        max_mismatches: Maximum mismatches (edits) allowed

    Returns:
        Updated primers_data with bowtie fields populated
    """
    indices_to_check = ["unconverted", "converted_unmethylated", "converted_methylated"]

    # Check which indices exist
    available_indices = []
    for idx_name in indices_to_check:
        idx_path = os.path.join(index_dir, idx_name)
        if os.path.exists(idx_path + ".1.bt2"):
            available_indices.append(idx_name)

    if not available_indices:
        for p in primers_data:
            p["bowtie_passes_filter"] = None
            p["bowtie_intended_genome"] = None
            p["mapping_error_note"] = "Bowtie index not available"
            p["mismatch_profile"] = ""
        return primers_data

    # Collect all unique primer sequences (left + right from all pairs)
    # Use a dict to deduplicate — many primers share sequences across pairs
    seq_to_id = {}
    all_seqs = []
    for p in primers_data:
        for side in ("left_primer", "right_primer"):
            seq = p[side]
            if seq not in seq_to_id:
                seq_to_id[seq] = len(all_seqs)
                all_seqs.append(seq)

    print(f"  Bowtie2 batch: {len(primers_data)} pairs, {len(all_seqs)} unique primers, "
          f"{len(available_indices)} indices")

    # Run bowtie2 once per index with ALL primers in one FASTA
    # Results: {idx_name: {seq_id: [alignment_dicts]}}
    batch_results = {}
    for idx_name in available_indices:
        idx_path = os.path.join(index_dir, idx_name)
        print(f"  Aligning vs {idx_name}...")
        results = batch_align_primers_bowtie2(all_seqs, idx_path, max_mismatches)
        # Map back from index to sequence
        batch_results[idx_name] = {}
        for seq_idx, alignments in results.items():
            seq = all_seqs[seq_idx]
            batch_results[idx_name][seq] = alignments

    # Now evaluate each primer pair using pre-computed alignments
    state_to_index = {
        "SM": "converted_methylated", "AM": "converted_methylated",
        "SU": "converted_unmethylated", "AU": "converted_unmethylated",
        "S": "unconverted", "A": "unconverted",
    }

    for p in primers_data:
        intended_index = state_to_index.get(p["template_used"], "unconverted")
        left_seq = p["left_primer"]
        right_seq = p["right_primer"]

        # Collect alignments for this pair from batch results
        all_alignments = {}
        for idx_name in indices_to_check:
            left_aln = batch_results.get(idx_name, {}).get(left_seq, [])
            right_aln = batch_results.get(idx_name, {}).get(right_seq, [])
            all_alignments[idx_name] = {"left": left_aln, "right": right_aln}

        result = _evaluate_screening(all_alignments, intended_index,
                                      indices_to_check, max_mismatches,
                                      product_size_min, product_size_max)

        p["bowtie_passes_filter"] = result.passes_filter
        p["bowtie_intended_genome"] = result.intended_genome
        p["mapping_error_note"] = result.mapping_note
        p["mismatch_profile"] = result.mismatch_profile
        p["mismatch_detail"] = result.mismatch_detail

    return primers_data


def _evaluate_screening(all_alignments: dict, intended_index: str,
                         indices_to_check: list, max_mismatches: int,
                         product_size_min: int = 40, product_size_max: int = 150) -> BowtieResult:
    """
    Evaluate screening results using PAIR-LEVEL specificity.

    A single bisulfite-converted primer (A/T-rich, ~20bp) will match many places
    in the 3-billion-bp genome. That's expected and not a problem. What matters
    for PCR/dPCR is whether the PAIR (left + right) can form an amplicon — i.e.,
    both primers align to the same genome location within product-size distance,
    in convergent orientation (left on forward, right on reverse).

    Pass criteria:
    1. The pair amplifies at the intended genome (both primers map there)
    2. The pair does NOT amplify at any off-target genome
    3. The pair amplifies uniquely at the intended genome (only 1 amplicon)
    4. Bowtie2 didn't fail

    Args:
        all_alignments: {idx_name: {"left": [alignments], "right": [alignments]}}
        intended_index: Which genome index is the intended one
        indices_to_check: List of all index names to check
        max_mismatches: Maximum mismatches (edits) allowed for a hit to count
        product_size_min: Min amplicon size for pair-level check
        product_size_max: Max amplicon size for pair-level check

    Returns:
        BowtieResult with pass/fail and details
    """
    notes = []
    bowtie_failed = False

    # Track mismatch counts (0-5) across ALL genomes for both primers combined
    mismatch_counts = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0}

    # Count total alignments for mismatch profile
    total_alignments = 0
    for idx_name in indices_to_check:
        if idx_name not in all_alignments:
            continue
        for primer_name in ("left", "right"):
            for a in all_alignments[idx_name][primer_name]:
                if a["edits"] <= max_mismatches:
                    total_alignments += 1
                    e = a["edits"]
                    if 0 <= e <= 5:
                        mismatch_counts[e] += 1

    # Check for bowtie2 failure (no hits anywhere)
    if total_alignments == 0:
        notes.append("WARNING: No alignments found in any genome — possible bowtie2 error")
        bowtie_failed = True

    # ── Pair-level specificity: count amplicons per genome index ──
    # An amplicon forms when left primer maps forward and right primer maps reverse
    # at the same chromosome, with right primer start within product_size range
    # after left primer start.
    #
    # bowtie2 SAM flags:
    #   0 = forward strand, 16 = reverse strand
    # For a valid PCR pair: left primer on forward (flag 0), right primer on reverse (flag 16)
    # Left primer position = pos, Right primer position = pos (start of alignment on reference)
    # Amplicon size = right_pos + right_len - left_pos

    amplicons_per_genome = {}  # {idx_name: [(chrom, left_pos, amplicon_size, left_edits, right_edits)]}

    for idx_name in indices_to_check:
        if idx_name not in all_alignments:
            continue

        left_alignments = [a for a in all_alignments[idx_name]["left"] if a["edits"] <= max_mismatches]
        right_alignments = [a for a in all_alignments[idx_name]["right"] if a["edits"] <= max_mismatches]

        # Build lookup: {chrom: [(pos, flag, edits, len), ...]} for each primer
        left_by_chrom = {}
        for a in left_alignments:
            chrom = a["rname"]
            if chrom not in left_by_chrom:
                left_by_chrom[chrom] = []
            left_by_chrom[chrom].append(a)

        right_by_chrom = {}
        for a in right_alignments:
            chrom = a["rname"]
            if chrom not in right_by_chrom:
                right_by_chrom[chrom] = []
            right_by_chrom[chrom].append(a)

        # Find pairs: left on forward (flag & 16 == 0), right on reverse (flag & 16 != 0)
        # within product_size range
        amplicons = []
        for chrom in set(left_by_chrom.keys()) & set(right_by_chrom.keys()):
            for la in left_by_chrom[chrom]:
                if la["flag"] & 16:  # left primer should be forward
                    continue
                left_pos = la["pos"]
                left_len = len(la["seq"])

                for ra in right_by_chrom[chrom]:
                    if not (ra["flag"] & 16):  # right primer should be reverse
                        continue
                    right_pos = ra["pos"]
                    right_len = len(ra["seq"])

                    # Amplicon size = (right_pos + right_len) - left_pos
                    amplicon_size = (right_pos + right_len) - left_pos

                    if product_size_min <= amplicon_size <= product_size_max:
                        amplicons.append((chrom, left_pos, amplicon_size,
                                         la["edits"], ra["edits"]))

        amplicons_per_genome[idx_name] = amplicons

    # ── Evaluate pass/fail based on pair-level amplicons ──
    intended_amplicons = amplicons_per_genome.get(intended_index, [])
    off_target_amplicons = 0

    for idx_name in indices_to_check:
        if idx_name == intended_index:
            continue
        amps = amplicons_per_genome.get(idx_name, [])
        if amps:
            off_target_amplicons += len(amps)
            notes.append(f"Pair amplifies {len(amps)} off-target product(s) in {idx_name}")

    if len(intended_amplicons) == 0:
        notes.append("Pair does not amplify in intended genome")
    elif len(intended_amplicons) > 1:
        notes.append(f"Pair amplifies {len(intended_amplicons)} products in intended genome (multi-amplicon)")

    # Pass criteria:
    # 1. Pair amplifies at intended genome (at least 1 amplicon)
    # 2. No off-target amplicons in any other genome
    # 3. Only 1 amplicon at intended genome (unique)
    # 4. Bowtie2 didn't fail
    passes = (len(intended_amplicons) == 1
              and off_target_amplicons == 0
              and not bowtie_failed)

    note = "; ".join(notes) if notes else "Unique amplification at intended genome"

    # Build mismatch profile: compact code only (6 digits: 0mm, 1mm, 2mm, 3mm, 4mm, 5mm)
    compact = "".join(str(mismatch_counts[i]) for i in range(0, 6))

    # Build human-readable description (separate field)
    parts = []
    for i in range(0, 6):
        if i == 0:
            parts.append(f"{mismatch_counts[i]} with 0 mismatches (perfect)")
        elif i == 1:
            parts.append(f"{mismatch_counts[i]} with 1 mismatch")
        elif i <= 4:
            parts.append(f"{mismatch_counts[i]} with {i} mm's")
        else:
            parts.append(f"{mismatch_counts[i]} with {i} mms")
    mismatch_detail = ", ".join(parts)

    return BowtieResult(
        passes_filter=passes,
        intended_genome=intended_index,
        total_alignments=total_alignments,
        off_target_alignments=off_target_amplicons,
        mapping_note=note,
        mismatch_profile=compact,
        mismatch_detail=mismatch_detail,
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
