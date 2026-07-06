#!/usr/bin/env python3
"""
Phase 1: Load DMR blocks and per-CpG data from the Excel output.

Reads the DMR_percpg_full_atlas_all_cell_types.xlsx file and prepares
DMR blocks with their per-CpG methylation data for primer design.

Each block is loaded with:
- Genomic coordinates (chr, start, end)
- CpG positions and labels
- Target and background methylation per CpG
- Cleanliness score
- Gene annotation
"""

import os
import re
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Optional, Dict


@dataclass
class CpGSite:
    """A single CpG site within a DMR block."""
    label: str           # CpG_1, CpG_2, ...
    position: int        # Genomic position (0-based)
    global_idx: int      # Global CpG index in the genome
    target_mean_beta: float
    background_mean_beta: float
    delta_beta: float


@dataclass
class DMRBlock:
    """A DMR block with per-CpG methylation data."""
    cell_type_id: str    # MONO, BCELL, NK, etc.
    seq_id: str          # MONO_0001, etc.
    rank: int
    chrom: str
    start: int           # 0-based start
    end: int             # 1-based end (BED format)
    num_cpgs: int
    block_len: int
    gene: str
    annotation: str
    direction: str
    target_mean_meth: float
    background_mean_meth: float
    delta_means: float
    delta_quants: float
    delta_maxmin: float
    ttest_pval: str
    mwtest_pval: str
    mvalue_ttest: str
    cleanliness_score: float
    target_celltype: str
    cpg_sites: List[CpGSite] = field(default_factory=list)

    @property
    def coordinates(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}"

    @property
    def cpg_positions(self) -> List[int]:
        return [c.position for c in self.cpg_sites]


def load_dmr_blocks(xlsx_path: str, cell_type_id: Optional[str] = None) -> List[DMRBlock]:
    """
    Load DMR blocks from the per-CpG Excel file.

    Args:
        xlsx_path: Path to DMR_percpg_full_atlas_all_cell_types.xlsx
        cell_type_id: If specified, only load blocks for this cell type

    Returns:
        List of DMRBlock objects with per-CpG data
    """
    # Read block summary
    block_df = pd.read_excel(xlsx_path, sheet_name="Block_Summary")
    if cell_type_id:
        block_df = block_df[block_df["cell_type_id"] == cell_type_id]

    # Read per-CpG data
    percpg_df = pd.read_excel(xlsx_path, sheet_name="Block_Summary")
    # Actually, per-CpG data is in per-cell-type sheets
    # We need to read from PerCpG_{ct} sheets

    blocks = []
    for _, row in block_df.iterrows():
        ct_id = row["cell_type_id"]
        seq_id = row["seq_id"]

        block = DMRBlock(
            cell_type_id=ct_id,
            seq_id=seq_id,
            rank=int(row["rank"]),
            chrom=row["Chr"],
            start=int(row["Start"]),
            end=int(row["End"]),
            num_cpgs=int(row["NumCpGs"]),
            block_len=int(row["BlockLen_bp"]),
            gene=str(row.get("Gene", "")),
            annotation=str(row.get("Annotation", "")),
            direction=str(row.get("Direction", "")),
            target_mean_meth=float(row["Target_Mean_Meth"]),
            background_mean_meth=float(row["Background_Mean_Meth"]),
            delta_means=float(row["Delta_Means"]),
            delta_quants=float(row["Delta_Quants"]),
            delta_maxmin=float(row["Delta_MaxMin"]),
            ttest_pval=str(row.get("Ttest_Pval", "")),
            mwtest_pval=str(row.get("MWtest_Pval", "")),
            mvalue_ttest=str(row.get("Mvalue_Ttest", "")),
            cleanliness_score=float(row.get("cleanliness_score", 0)),
            target_celltype=str(row.get("Target_CellType", "")),
        )
        blocks.append(block)

    # Load per-CpG data from PerCpG sheets
    sheet_name = f"PerCpG_{cell_type_id}" if cell_type_id else None

    if cell_type_id:
        # Read specific cell type sheet
        try:
            percpg_df = pd.read_excel(xlsx_path, sheet_name=f"PerCpG_{cell_type_id}")
        except Exception:
            percpg_df = None
    else:
        # Read all PerCpG sheets and concatenate
        all_percpg = []
        xl = pd.ExcelFile(xlsx_path)
        for sheet in xl.sheet_names:
            if sheet.startswith("PerCpG_"):
                df = pd.read_excel(xlsx_path, sheet_name=sheet)
                all_percpg.append(df)
        percpg_df = pd.concat(all_percpg, ignore_index=True) if all_percpg else None

    # Attach CpG sites to blocks
    if percpg_df is not None:
        for block in blocks:
            block_cpgs = percpg_df[percpg_df["seq_id"] == block.seq_id]
            for _, crow in block_cpgs.iterrows():
                cpg = CpGSite(
                    label=str(crow["cpg_label"]),
                    position=int(crow["cpg_position"]),
                    global_idx=int(crow["cpg_global_idx"]),
                    target_mean_beta=float(crow.get("target_mean_beta", 0) or 0),
                    background_mean_beta=float(crow.get("background_mean_beta", 0) or 0),
                    delta_beta=float(crow.get("delta_beta", 0) or 0),
                )
                block.cpg_sites.append(cpg)

    return blocks


def load_cff_regions(xlsx_path: str) -> pd.DataFrame:
    """Load CFF regions from the CFF Excel file."""
    sheets = pd.read_excel(xlsx_path, sheet_name=None)
    results = []
    for sheet_name, df in sheets.items():
        if sheet_name.startswith("CFF_"):
            df["assay_type"] = sheet_name
            results.append(df)
    if results:
        return pd.concat(results, ignore_index=True)
    return pd.DataFrame()


def load_conversion_controls(xlsx_path: str) -> pd.DataFrame:
    """Load conversion control candidates."""
    return pd.read_excel(xlsx_path, sheet_name="Conversion_Controls")


def filter_blocks_by_cpg_count(blocks: List[DMRBlock], min_cpg: int = 5,
                                preferred_min_cpg: int = 7) -> List[DMRBlock]:
    """Filter blocks by minimum CpG count, preferring blocks with more CpGs."""
    filtered = [b for b in blocks if b.num_cpgs >= min_cpg]
    # Sort by CpG count descending (prefer more CpGs)
    filtered.sort(key=lambda b: (-b.num_cpgs, -b.cleanliness_score))
    return filtered


def filter_blocks_by_cleanliness(blocks: List[DMRBlock], min_score: float = 0.6) -> List[DMRBlock]:
    """Filter blocks by cleanliness score."""
    return [b for b in blocks if b.cleanliness_score >= min_score]


def main():
    """CLI entry point for testing."""
    import argparse
    parser = argparse.ArgumentParser(description="Load DMR blocks from Excel")
    parser.add_argument("--xlsx", required=True, help="Path to DMR Excel file")
    parser.add_argument("--cell-type", help="Filter by cell type ID (e.g. MONO)")
    parser.add_argument("--min-cpg", type=int, default=5, help="Minimum CpGs per block")
    parser.add_argument("--min-score", type=float, default=0.6, help="Minimum cleanliness score")
    args = parser.parse_args()

    blocks = load_dmr_blocks(args.xlsx, args.cell_type)
    blocks = filter_blocks_by_cpg_count(blocks, args.min_cpg)
    blocks = filter_blocks_by_cleanliness(blocks, args.min_score)

    print(f"Loaded {len(blocks)} blocks")
    for b in blocks[:5]:
        print(f"  {b.seq_id}: {b.coordinates} {b.num_cpgs}CpGs score={b.cleanliness_score:.3f} gene={b.gene}")
        for c in b.cpg_sites[:3]:
            print(f"    {c.label} pos={c.position} target={c.target_mean_beta:.3f} bg={c.background_mean_beta:.3f}")


if __name__ == "__main__":
    main()
