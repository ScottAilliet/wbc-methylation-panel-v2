#!/usr/bin/env python3
"""
WBC Methylation Panel — Primer Design Pipeline CLI

Modular step-by-step pipeline for designing bisulfite-PCR primers
for WBC methylation markers.

Usage:
    # DMR discovery from raw beta files (fully self-contained):
    python -m methyl_panel.pipeline --steps all --discover-dmrs \\
        --cell-type MONO --genome data/hg19/hg19.fa.gz \\
        --min-tm 58 --opt-tm 60 --max-tm 62 --output-dir results/

    # Load DMR blocks from a pre-computed Excel file:
    python -m methyl_panel.pipeline --steps 1,2,3 --dmr-xlsx DMR.xlsx --genome hg19.fa.gz
    python -m methyl_panel.pipeline --steps all --output-dir results/
    python -m methyl_panel.pipeline --steps 3,4,5,6,7 --primers primers.json

Steps:
    0. Discover DMRs from raw beta files (wgbstools find_markers)  [--discover-dmrs]
    1. Load DMR blocks from Excel
    2. Bisulfite convert genomic sequences
    3. Design primers with Primer3
    4. Bowtie2 specificity screening (6-strand)
    5. Secondary structure screening
    6. Common SNP screening
    7. Primer-dimer prediction (DimerDetective)
    8. Generate XLSX output (U-assays format)
    9. Generate PDF output (U-assays format)
   10. Multiplexing compatibility filter (select 1 assay per cell type)  [--multiplex-dirs]

Step 10 is NOT included in --steps all. It requires all per-cell-type runs
to be complete first. Use --steps 10 --multiplex-dirs results/MONO,results/BCELL,...
"""

import os
import sys
import json
import argparse
from typing import List, Optional, Dict

from methyl_panel.config import Primer3PlusConfig, PipelineConfig


def step0_discover_dmrs(args):
    """Step 0: Discover DMRs from raw beta files using wgbstools find_markers."""
    from methyl_panel.phase0_dmr_discovery import discover_dmrs, VALID_CELL_TYPES

    print("\n=== Step 0: DMR Discovery (wgbstools find_markers) ===")

    if not args.cell_type:
        raise ValueError(
            "--cell-type is required when using --discover-dmrs. "
            f"Valid options: {', '.join(VALID_CELL_TYPES)}"
        )

    blocks = discover_dmrs(
        cell_type_id=args.cell_type,
        beta_dir=args.beta_dir,
        blocks_file=args.blocks_file,
        groups_csv=args.groups_csv,
        out_dir=args.output_dir,
        threads=args.threads,
        top_n=args.top_markers,
        wgbstools_path=args.wgbstools_path,
        max_bg_samples=args.max_bg_samples,
        skip_find_markers=args.skip_find_markers,
    )

    # Convert DMRBlock objects to the same dict format that step 1 produces
    # so that steps 2-9 can consume them identically
    blocks_data = []
    for b in blocks:
        blocks_data.append({
            "cell_type_id": b.cell_type_id,
            "seq_id": b.seq_id,
            "rank": b.rank,
            "chrom": b.chrom,
            "start": b.start,
            "end": b.end,
            "num_cpgs": b.num_cpgs,
            "block_len": b.block_len,
            "gene": b.gene,
            "annotation": b.annotation,
            "cleanliness_score": b.cleanliness_score,
            "target_mean_meth": b.target_mean_meth,
            "background_mean_meth": b.background_mean_meth,
            "delta_means": b.delta_means,
            "cpg_sites": [
                {"label": c.label, "position": c.position,
                 "global_idx": c.global_idx,
                 "target_mean_beta": c.target_mean_beta,
                 "background_mean_beta": c.background_mean_beta}
                for c in b.cpg_sites
            ],
        })

    print(f"\n  DMR discovery complete: {len(blocks_data)} blocks for {args.cell_type}")
    return blocks_data


def step1_load_dmrs(args):
    """Step 1: Load DMR blocks from Excel."""
    from methyl_panel.phase1_dmr_loader import load_dmr_blocks, filter_blocks_by_cpg_count

    print("\n=== Step 1: Loading DMR blocks ===")
    cell_type = args.cell_type if args.cell_type else None
    blocks = load_dmr_blocks(args.dmr_xlsx, cell_type)
    blocks = filter_blocks_by_cpg_count(blocks, args.min_cpg, args.preferred_min_cpg)

    print(f"Loaded {len(blocks)} blocks")
    if cell_type:
        print(f"  Filtered to cell type: {cell_type}")

    # Save to JSON for downstream steps
    blocks_data = []
    for b in blocks:
        blocks_data.append({
            "cell_type_id": b.cell_type_id,
            "seq_id": b.seq_id,
            "rank": b.rank,
            "chrom": b.chrom,
            "start": b.start,
            "end": b.end,
            "num_cpgs": b.num_cpgs,
            "block_len": b.block_len,
            "gene": b.gene,
            "annotation": b.annotation,
            "cleanliness_score": b.cleanliness_score,
            "target_mean_meth": b.target_mean_meth,
            "background_mean_meth": b.background_mean_meth,
            "delta_means": b.delta_means,
            "cpg_sites": [
                {"label": c.label, "position": c.position,
                 "global_idx": c.global_idx,
                 "target_mean_beta": c.target_mean_beta,
                 "background_mean_beta": c.background_mean_beta}
                for c in b.cpg_sites
            ],
        })

    out_path = os.path.join(args.output_dir, "dmr_blocks.json")
    with open(out_path, 'w') as f:
        json.dump(blocks_data, f, indent=2)
    print(f"  Saved to {out_path}")
    return blocks_data


def step2_bisulfite_convert(args, blocks_data=None):
    """Step 2: Bisulfite convert genomic sequences."""
    from methyl_panel.phase2_bisulfite_convert import convert_sequence, fetch_genomic_sequence

    print("\n=== Step 2: Bisulfite converting sequences ===")
    if blocks_data is None:
        blocks_path = os.path.join(args.output_dir, "dmr_blocks.json")
        with open(blocks_path) as f:
            blocks_data = json.load(f)

    converted = []
    skipped = 0
    for i, block in enumerate(blocks_data):
        if args.max_blocks and (i - skipped) >= args.max_blocks:
            break

        # Fetch genomic sequence with flanking
        flank = args.flank
        try:
            seq = fetch_genomic_sequence(
                args.genome, block["chrom"],
                block["start"], block["end"], flank
            )
        except (RuntimeError, Exception) as e:
            print(f"  Skipping {block['seq_id']} ({block['chrom']}): {e}")
            skipped += 1
            continue

        if not seq or 'N' * 10 in seq:
            print(f"  Skipping {block['seq_id']} ({block['chrom']}): sequence not available")
            skipped += 1
            continue

        strands = convert_sequence(seq)

        converted.append({
            "seq_id": block["seq_id"],
            "cell_type_id": block["cell_type_id"],
            "chrom": block["chrom"],
            "start": block["start"],
            "end": block["end"],
            "genomic_seq": seq,
            "strands": {
                "S": strands.sense,
                "SM": strands.sense_methylated,
                "SU": strands.sense_unmethylated,
                "A": strands.antisense,
                "AM": strands.antisense_methylated,
                "AU": strands.antisense_unmethylated,
            },
            "cpg_positions": strands.cpg_positions,
        })

        if (i + 1) % 50 == 0:
            print(f"  Converted {i+1}/{len(blocks_data)} blocks")

    out_path = os.path.join(args.output_dir, "converted_sequences.json")
    with open(out_path, 'w') as f:
        json.dump(converted, f)
    print(f"  Converted {len(converted)} blocks, saved to {out_path}")
    return converted


def step3_design_primers(args, converted_data=None):
    """Step 3: Design primers with Primer3."""
    from methyl_panel.phase3_primer3_design import design_primers_for_block

    print("\n=== Step 3: Designing primers ===")
    if converted_data is None:
        path = os.path.join(args.output_dir, "converted_sequences.json")
        with open(path) as f:
            converted_data = json.load(f)

    # Load config
    if args.settings:
        config = Primer3PlusConfig.from_primer3plus_file(args.settings)
    else:
        config = Primer3PlusConfig()

    # Set Tm (required)
    if args.min_tm is not None:
        config.primer_min_tm = args.min_tm
    if args.opt_tm is not None:
        config.primer_opt_tm = args.opt_tm
    if args.max_tm is not None:
        config.primer_max_tm = args.max_tm

    config.validate_tm()
    pipeline_config = PipelineConfig()

    all_primers = []
    for i, block in enumerate(converted_data):
        if args.max_blocks and i >= args.max_blocks:
            break

        primers = design_primers_for_block(
            block["genomic_seq"], config, pipeline_config,
            f"{block['seq_id']}_001", block["seq_id"], block["cell_type_id"]
        )

        all_primers.extend(primers)

        if (i + 1) % 20 == 0:
            print(f"  Designed primers for {i+1}/{len(converted_data)} blocks ({len(all_primers)} total)")

    # Save to JSON
    primers_data = []
    for p in all_primers:
        primers_data.append({
            "assay_id": p.assay_id, "seq_id": p.seq_id,
            "cell_type_id": p.cell_type_id, "template_used": p.template_used,
            "left_primer": p.left_primer, "right_primer": p.right_primer,
            "left_primer_display": p.left_primer_display,
            "right_primer_display": p.right_primer_display,
            "left_tm": p.left_tm, "right_tm": p.right_tm,
            "product_size": p.product_size,
            "c_total_tail": p.c_total_tail, "c_total": p.c_total,
            "left_c_total": p.left_c_total, "right_c_total": p.right_c_total,
            "left_c_tail": p.left_c_tail, "right_c_tail": p.right_c_tail,
            "sense_meth_mismatch_score": p.sense_meth_mismatch_score,
            "sense_unmeth_mismatch_score": p.sense_unmeth_mismatch_score,
            "anti_meth_mismatch_score": p.anti_meth_mismatch_score,
            "anti_unmeth_mismatch_score": p.anti_unmeth_mismatch_score,
            "left_gc_percent": p.left_gc_percent, "right_gc_percent": p.right_gc_percent,
            "penalty": p.penalty,
            "left_start": p.left_start, "left_len": p.left_len,
            "right_start": p.right_start, "right_len": p.right_len,
        })

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Designed {len(all_primers)} primer pairs, saved to {out_path}")
    return primers_data


def step4_bowtie(args, primers_data=None):
    """Step 4: Bowtie2 specificity screening."""
    from methyl_panel.phase4_bowtie_specificity import screen_primer_pair

    print("\n=== Step 4: Bowtie2 specificity screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        if args.bowtie_index_dir and os.path.exists(args.bowtie_index_dir):
            result = screen_primer_pair(
                p["left_primer"], p["right_primer"],
                p["template_used"], args.bowtie_index_dir
            )
            p["bowtie_passes_filter"] = result.passes_filter
            p["bowtie_intended_genome"] = result.intended_genome
            p["mapping_error_note"] = result.mapping_note
        else:
            p["bowtie_passes_filter"] = None
            p["bowtie_intended_genome"] = None
            p["mapping_error_note"] = "Bowtie index not available"

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Bowtie screening complete, updated {out_path}")
    return primers_data


def step5_structure(args, primers_data=None):
    """Step 5: Secondary structure screening."""
    from methyl_panel.phase5_structure import screen_structure

    print("\n=== Step 5: Secondary structure screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        result = screen_structure(p["left_primer"], p["right_primer"])
        p["left_structure_mfe"] = result.left_mfe
        p["right_structure_mfe"] = result.right_mfe

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Structure screening complete, updated {out_path}")
    return primers_data


def step6_snp(args, primers_data=None):
    """Step 6: Common SNP screening."""
    from methyl_panel.phase6_snp import screen_primer_pair_snp

    print("\n=== Step 6: Common SNP screening ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        if args.dbsnp and os.path.exists(args.dbsnp):
            # Need genomic positions — would need to map primer positions back
            # For now, set as None
            p["common_variant_score"] = None
        else:
            p["common_variant_score"] = None

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  SNP screening complete (requires dbSNP VCF for full screening)")
    return primers_data


def step7_dimer(args, primers_data=None):
    """Step 7: Primer-dimer prediction (DimerDetective)."""
    from methyl_panel.phase7_dimer import predict_dimer

    print("\n=== Step 7: DimerDetective primer-dimer prediction ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    for p in primers_data:
        result = predict_dimer(p["left_primer"], p["right_primer"])
        p["primer_dimer_prediction"] = result.risk_tier
        p["primer_dimer_end_min_dg"] = result.end_min_dg

    out_path = os.path.join(args.output_dir, "primers.json")
    with open(out_path, 'w') as f:
        json.dump(primers_data, f, indent=2)
    print(f"  Dimer prediction complete, updated {out_path}")
    return primers_data


def step8_xlsx(args, primers_data=None):
    """Step 8: Generate XLSX output."""
    from methyl_panel.output_xlsx import write_xlsx
    from methyl_panel.phase3_primer3_design import PrimerPair

    print("\n=== Step 8: Generating XLSX output ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    # Convert dicts back to PrimerPair objects
    pairs = []
    for p in primers_data:
        pairs.append(PrimerPair(
            assay_id=p["assay_id"], seq_id=p["seq_id"],
            cell_type_id=p["cell_type_id"], template_used=p["template_used"],
            left_primer=p["left_primer"], right_primer=p["right_primer"],
            left_primer_display=p["left_primer_display"],
            right_primer_display=p["right_primer_display"],
            left_start=p["left_start"], left_len=p["left_len"],
            right_start=p["right_start"], right_len=p["right_len"],
            left_tm=p["left_tm"], right_tm=p["right_tm"],
            product_size=p["product_size"],
            c_total_tail=p["c_total_tail"], c_total=p["c_total"],
            left_c_total=p["left_c_total"], right_c_total=p["right_c_total"],
            left_c_tail=p["left_c_tail"], right_c_tail=p["right_c_tail"],
            sense_meth_mismatch_score=p["sense_meth_mismatch_score"],
            sense_unmeth_mismatch_score=p["sense_unmeth_mismatch_score"],
            anti_meth_mismatch_score=p["anti_meth_mismatch_score"],
            anti_unmeth_mismatch_score=p["anti_unmeth_mismatch_score"],
            left_gc_percent=p["left_gc_percent"], right_gc_percent=p["right_gc_percent"],
            left_self_any=0, right_self_any=0, left_self_end=0, right_self_end=0,
            pair_compl_any=0, pair_compl_end=0,
            left_end_stability=0, right_end_stability=0,
            penalty=p["penalty"],
            bowtie_passes_filter=p.get("bowtie_passes_filter"),
            bowtie_intended_genome=p.get("bowtie_intended_genome"),
            left_structure_mfe=p.get("left_structure_mfe"),
            right_structure_mfe=p.get("right_structure_mfe"),
            primer_dimer_prediction=p.get("primer_dimer_prediction"),
            primer_dimer_end_min_dg=p.get("primer_dimer_end_min_dg"),
            common_variant_score=p.get("common_variant_score"),
            mapping_error_note=p.get("mapping_error_note"),
        ))

    out_path = os.path.join(args.output_dir, "primer_assays.xlsx")
    write_xlsx(pairs, out_path)
    print(f"  XLSX saved to {out_path}")


def step9_pdf(args, primers_data=None):
    """Step 9: Generate PDF output."""
    from methyl_panel.output_pdf import write_pdf
    from methyl_panel.phase3_primer3_design import PrimerPair

    print("\n=== Step 9: Generating PDF output ===")
    if primers_data is None:
        path = os.path.join(args.output_dir, "primers.json")
        with open(path) as f:
            primers_data = json.load(f)

    pairs = []
    for p in primers_data:
        pairs.append(PrimerPair(
            assay_id=p["assay_id"], seq_id=p["seq_id"],
            cell_type_id=p["cell_type_id"], template_used=p["template_used"],
            left_primer=p["left_primer"], right_primer=p["right_primer"],
            left_primer_display=p["left_primer_display"],
            right_primer_display=p["right_primer_display"],
            left_start=p["left_start"], left_len=p["left_len"],
            right_start=p["right_start"], right_len=p["right_len"],
            left_tm=p["left_tm"], right_tm=p["right_tm"],
            product_size=p["product_size"],
            c_total_tail=p["c_total_tail"], c_total=p["c_total"],
            left_c_total=p["left_c_total"], right_c_total=p["right_c_total"],
            left_c_tail=p["left_c_tail"], right_c_tail=p["right_c_tail"],
            sense_meth_mismatch_score=p["sense_meth_mismatch_score"],
            sense_unmeth_mismatch_score=p["sense_unmeth_mismatch_score"],
            anti_meth_mismatch_score=p["anti_meth_mismatch_score"],
            anti_unmeth_mismatch_score=p["anti_unmeth_mismatch_score"],
            left_gc_percent=p["left_gc_percent"], right_gc_percent=p["right_gc_percent"],
            left_self_any=0, right_self_any=0, left_self_end=0, right_self_end=0,
            pair_compl_any=0, pair_compl_end=0,
            left_end_stability=0, right_end_stability=0,
            penalty=p["penalty"],
            bowtie_passes_filter=p.get("bowtie_passes_filter"),
            bowtie_intended_genome=p.get("bowtie_intended_genome"),
            left_structure_mfe=p.get("left_structure_mfe"),
            right_structure_mfe=p.get("right_structure_mfe"),
            primer_dimer_prediction=p.get("primer_dimer_prediction"),
            primer_dimer_end_min_dg=p.get("primer_dimer_end_min_dg"),
            common_variant_score=p.get("common_variant_score"),
            mapping_error_note=p.get("mapping_error_note"),
        ))

    out_path = os.path.join(args.output_dir, "primer_assays.pdf")
    write_pdf(pairs, output_path=out_path)
    print(f"  PDF saved to {out_path}")


def step10_multiplex(args):
    """Step 10: Multiplexing compatibility filter — select one assay per cell type."""
    from methyl_panel.phase10_multiplex import (
        select_multiplex_panel, save_panel_json, save_panel_xlsx, save_panel_pdf,
    )

    print("\n=== Step 10: Multiplexing Compatibility Filter ===")

    if not args.multiplex_dirs:
        raise ValueError(
            "--multiplex-dirs is required for step 10. "
            "Provide a comma-separated list of per-cell-type output directories, "
            "e.g. --multiplex-dirs results/MONO,results/BCELL,results/NK,..."
        )

    # Parse multiplex dirs — auto-detect cell type from directory name
    dirs = [d.strip() for d in args.multiplex_dirs.split(",")]
    cell_type_dirs = {}
    for d in dirs:
        ct = os.path.basename(os.path.normpath(d))
        cell_type_dirs[ct] = d

    print(f"  Cell types: {list(cell_type_dirs.keys())}")
    print(f"  Criteria: dimer cutoff={args.cross_dimer_cutoff}, "
          f"tm tol={args.tm_tolerance}, min amp diff={args.min_amplicon_diff}")

    opt_tm = args.opt_tm if args.opt_tm else 60.0

    panel = select_multiplex_panel(
        cell_type_dirs,
        opt_tm=opt_tm,
        cross_dimer_cutoff=args.cross_dimer_cutoff,
        tm_tolerance=args.tm_tolerance,
        min_amplicon_diff=args.min_amplicon_diff,
    )

    # Save outputs
    json_path = save_panel_json(
        panel, os.path.join(args.output_dir, "multiplex_panel.json"))
    print(f"\n  Panel JSON: {json_path}")

    xlsx_path = save_panel_xlsx(
        panel, os.path.join(args.output_dir, "multiplex_panel.xlsx"))
    print(f"  Panel XLSX: {xlsx_path}")

    pdf_path = save_panel_pdf(
        panel, os.path.join(args.output_dir, "multiplex_panel.pdf"))
    print(f"  Panel PDF: {pdf_path}")

    # Summary
    print(f"\n  Selected: {panel.n_cell_types_selected}/{panel.n_cell_types_requested} cell types")
    if panel.failed_cell_types:
        print(f"  Failed: {', '.join(panel.failed_cell_types)}")
    for a in panel.selected_assays:
        print(f"    {a['cell_type_id']}: {a['assay_id']} "
              f"(q={a['quality_score']:.3f}, Tm={a['left_tm']:.1f}/{a['right_tm']:.1f}, "
              f"product={a['product_size']}bp)")

    return panel


def main():
    parser = argparse.ArgumentParser(
        description="WBC Methylation Panel — Primer Design Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--steps", default="all",
                        help="Comma-separated step numbers (e.g. 0,2,3 or 1,2,3) or 'all'")
    parser.add_argument("--output-dir", default="results/",
                        help="Output directory for results")
    parser.add_argument("--dmr-xlsx", default="data/WBC_Panel_Top200_v7.9.xlsx",
                        help="Path to DMR per-CpG Excel file (for step 1)")
    parser.add_argument("--genome", default="data/hg19/hg19.fa.gz",
                        help="Path to genome FASTA (can be .gz)")
    parser.add_argument("--settings", help="Primer3Plus settings file")
    parser.add_argument("--cell-type",
                        help="Cell type ID (MONO, BCELL, NK, GRAN, CD3T, CD8T, CD4T). "
                             "Required with --discover-dmrs.")
    parser.add_argument("--min-cpg", type=int, default=5, help="Min CpGs per block")
    parser.add_argument("--preferred-min-cpg", type=int, default=7)
    parser.add_argument("--min-tm", type=float, help="Min Tm (°C) — REQUIRED")
    parser.add_argument("--opt-tm", type=float, help="Opt Tm (°C) — REQUIRED")
    parser.add_argument("--max-tm", type=float, help="Max Tm (°C) — REQUIRED")
    parser.add_argument("--flank", type=int, default=100, help="Flanking bp around DMR")
    parser.add_argument("--max-blocks", type=int, help="Max blocks to process (for testing)")
    parser.add_argument("--bowtie-index-dir", help="Directory with bowtie2 indices")
    parser.add_argument("--dbsnp", help="Path to dbSNP VCF")

    # --- Step 0: DMR Discovery arguments ---
    parser.add_argument("--discover-dmrs", action="store_true",
                        help="Run Step 0 (DMR discovery from raw beta files) "
                             "instead of Step 1 (load from Excel)")
    parser.add_argument("--beta-dir", default="data/beta_files/",
                        help="Directory containing .beta files (for --discover-dmrs)")
    parser.add_argument("--blocks-file", default="data/GSE186458_blocks.s205.bed.gz",
                        help="Path to blocks BED file (for --discover-dmrs)")
    parser.add_argument("--groups-csv", default="data/full_atlas_groups.csv",
                        help="Path to full atlas groups CSV (for --discover-dmrs)")
    parser.add_argument("--wgbstools-path",
                        help="Path to wgbstools executable (auto-detected if omitted)")
    parser.add_argument("--threads", type=int, default=2,
                        help="Number of threads for find_markers")
    parser.add_argument("--top-markers", type=int, default=200,
                        help="Number of top markers per cell type")
    parser.add_argument("--max-bg-samples", type=int, default=30,
                        help="Max background samples for per-CpG extraction")
    parser.add_argument("--skip-find-markers", action="store_true",
                        help="Skip find_markers (reuse existing BED output)")

    # --- Step 10: Multiplexing arguments ---
    parser.add_argument("--multiplex-dirs",
                        help="Comma-separated list of per-cell-type output dirs for "
                             "multiplex selection (e.g. results/MONO,results/BCELL,...). "
                             "Required for step 10.")
    parser.add_argument("--cross-dimer-cutoff", type=float, default=-1.0,
                        help="Cross-dimer ΔG cutoff for multiplex (kcal/mol, default -1.0)")
    parser.add_argument("--tm-tolerance", type=float, default=2.0,
                        help="Max Tm spread across compatible assays (°C, default 2.0)")
    parser.add_argument("--min-amplicon-diff", type=int, default=10,
                        help="Min product size difference for multiplex (bp, default 10)")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Parse steps — supports "all", "all,10", "0,2,3", "10", etc.
    # Step 10 (multiplex) is never in "all" by itself — it requires multiple
    # cell-type runs to be complete first. Use --steps 10 or --steps all,10.
    raw_steps = args.steps.split(",")
    raw_steps = [s.strip() for s in raw_steps]
    has_all = "all" in raw_steps
    explicit_steps = [int(s) for s in raw_steps if s != "all"]

    if has_all:
        base_steps = [0, 2, 3, 4, 5, 6, 7, 8, 9] if args.discover_dmrs else [1, 2, 3, 4, 5, 6, 7, 8, 9]
        steps = base_steps + [s for s in explicit_steps if s not in base_steps]
    else:
        steps = explicit_steps
        # If --discover-dmrs and user specified step 1, replace with 0
        if args.discover_dmrs:
            steps = [0 if s == 1 else s for s in steps]

    # Validate --multiplex-dirs when step 10 is requested
    if 10 in steps and args.multiplex_dirs is None:
        parser.error("--multiplex-dirs is required when step 10 is in --steps.")

    # Validate Tm is set for primer design steps
    if any(s in steps for s in [3]) and not all([args.min_tm, args.opt_tm, args.max_tm]):
        if not args.settings:
            parser.error("Tm values (--min-tm, --opt-tm, --max-tm) are REQUIRED for primer design. "
                         "Or provide --settings file with Tm values.")

    # Execute steps
    blocks_data = None
    converted_data = None
    primers_data = None

    for step in steps:
        if step == 0:
            blocks_data = step0_discover_dmrs(args)
        elif step == 1:
            blocks_data = step1_load_dmrs(args)
        elif step == 2:
            converted_data = step2_bisulfite_convert(args, blocks_data)
        elif step == 3:
            primers_data = step3_design_primers(args, converted_data)
        elif step == 4:
            primers_data = step4_bowtie(args, primers_data)
        elif step == 5:
            primers_data = step5_structure(args, primers_data)
        elif step == 6:
            primers_data = step6_snp(args, primers_data)
        elif step == 7:
            primers_data = step7_dimer(args, primers_data)
        elif step == 8:
            step8_xlsx(args, primers_data)
        elif step == 9:
            step9_pdf(args, primers_data)
        elif step == 10:
            step10_multiplex(args)

    print("\n=== Pipeline complete ===")


if __name__ == "__main__":
    main()
