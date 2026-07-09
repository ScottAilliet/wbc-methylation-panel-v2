#!/usr/bin/env python3
"""
Phase 1: Load DMR blocks and per-CpG data from the Excel output.

Supports two Excel formats:

1. **WBC Panel v7.9 format** (one sheet per cell type):
   Sheets named CD4T, CD8T, B, NK, Mono, Gran, Blood-T-CD3.
   Each sheet has columns: Rank, Chr, Start, End, NumCpGs, BlockLen_bp,
   Gene, Annotation, Direction, Target_Mean_Meth, Background_Mean_Meth,
   Delta_Means, ..., Cleanliness_Score, tg_cpg1_mean, bg_cpg1_mean, ...
   Per-CpG positions are in block_cpg_coords / block_cpg_labels columns.

2. **Block_Summary format** (original v2 design):
   A Block_Summary sheet with cell_type_id, seq_id, Chr, Start, End, ...
   plus PerCpG_{cell_type_id} sheets with per-CpG data.

The loader auto-detects which format is present.
"""

import os
import re
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple


# ---------------------------------------------------------------------------
# Cell type ID mapping
# ---------------------------------------------------------------------------
# v7.9 sheet name -> short cell type ID used throughout the pipeline
SHEET_TO_CTID: Dict[str, str] = {
    "Mono": "MONO",
    "B": "BCELL",
    "NK": "NK",
    "Gran": "GRAN",
    "Blood-T-CD3": "CD3T",
    "CD8T": "CD8T",
    "CD4T": "CD4T",
}

# Reverse mapping: short ID -> v7.9 sheet name
CTID_TO_SHEET: Dict[str, str] = {v: k for k, v in SHEET_TO_CTID.items()}

# Sheets to skip (not cell-type data sheets)
NON_DATA_SHEETS = {"Glossary", "Summary", "Analysis_Summary"}


@dataclass
class CpGSite:
    """A single CpG site within a DMR block."""
    label: str           # CpG_1, CpG_2, ...
    position: int        # Genomic position (1-based)
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
    start: int           # 1-based start (inclusive)
    end: int             # 1-based end (inclusive)
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


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------
def _detect_format(xlsx_path: str) -> str:
    """Detect whether the Excel file is v7.9 or Block_Summary format."""
    xl = pd.ExcelFile(xlsx_path)
    sheet_names = set(xl.sheet_names)
    if "Block_Summary" in sheet_names:
        return "block_summary"
    # Check if any sheet matches v7.9 cell type names
    v79_sheets = set(SHEET_TO_CTID.keys()) | {"CD4T", "CD8T"}
    if sheet_names & v79_sheets:
        return "v79"
    # Fallback: check if first non-skip sheet has 'Rank' and 'Chr' columns
    for sn in xl.sheet_names:
        if sn in NON_DATA_SHEETS:
            continue
        df = pd.read_excel(xlsx_path, sheet_name=sn, nrows=1)
        if "Rank" in df.columns and "Chr" in df.columns:
            return "v79"
    return "unknown"


# ---------------------------------------------------------------------------
# v7.9 format loader
# ---------------------------------------------------------------------------
def _parse_cpg_coords(coords_str: str) -> List[Tuple[str, int]]:
    """Parse 'chr10:76391809, chr10:76391817' into [(chr, pos), ...]."""
    if not coords_str or pd.isna(coords_str):
        return []
    results = []
    for part in str(coords_str).split(","):
        part = part.strip()
        if ":" not in part:
            continue
        chrom, pos_str = part.rsplit(":", 1)
        try:
            pos = int(pos_str)
        except ValueError:
            continue
        results.append((chrom, pos))
    return results


def _load_v79(xlsx_path: str, cell_type_id: Optional[str] = None) -> List[DMRBlock]:
    """Load DMR blocks from a v7.9-format Excel file."""
    xl = pd.ExcelFile(xlsx_path)

    # Determine which sheets to read
    if cell_type_id:
        # Map short ID to sheet name
        sheet_name = CTID_TO_SHEET.get(cell_type_id, cell_type_id)
        if sheet_name not in xl.sheet_names:
            raise ValueError(
                f"Cell type '{cell_type_id}' not found. "
                f"Available sheets: {[s for s in xl.sheet_names if s not in NON_DATA_SHEETS]}"
            )
        sheets_to_read = [sheet_name]
    else:
        sheets_to_read = [s for s in xl.sheet_names if s not in NON_DATA_SHEETS]

    blocks = []
    for sheet_name in sheets_to_read:
        ct_id = SHEET_TO_CTID.get(sheet_name, sheet_name)
        df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

        for _, row in df.iterrows():
            rank = int(row["Rank"])
            seq_id = f"{ct_id}_{rank:04d}"

            chrom = str(row["Chr"])
            start = int(row["Start"])
            end = int(row["End"])
            num_cpgs = int(row["NumCpGs"])

            # Parse CpG positions from block_cpg_coords
            coords_str = row.get("block_cpg_coords", "")
            cpg_coords = _parse_cpg_coords(coords_str)

            # If block_cpg_coords is missing, distribute CpGs evenly across the block
            if not cpg_coords and num_cpgs > 0:
                step = max(1, (end - start) // max(1, num_cpgs))
                cpg_coords = [(chrom, start + i * step) for i in range(num_cpgs)]

            # Build CpG sites
            cpg_sites = []
            for i, (cpg_chrom, cpg_pos) in enumerate(cpg_coords, 1):
                # Get per-CpG methylation from tg_cpgN_mean / bg_cpgN_mean columns
                tg_col = f"tg_cpg{i}_mean"
                bg_col = f"bg_cpg{i}_mean"
                tg_val = float(row[tg_col]) if tg_col in row and pd.notna(row[tg_col]) else 0.0
                bg_val = float(row[bg_col]) if bg_col in row and pd.notna(row[bg_col]) else 0.0

                cpg_sites.append(CpGSite(
                    label=f"CpG_{i}",
                    position=cpg_pos,
                    global_idx=cpg_pos,  # Use genomic position as global index
                    target_mean_beta=tg_val,
                    background_mean_beta=bg_val,
                    delta_beta=abs(tg_val - bg_val),
                ))

            block = DMRBlock(
                cell_type_id=ct_id,
                seq_id=seq_id,
                rank=rank,
                chrom=chrom,
                start=start,
                end=end,
                num_cpgs=num_cpgs,
                block_len=int(row.get("BlockLen_bp", end - start + 1)),
                gene=str(row.get("Gene", "")),
                annotation=str(row.get("Annotation", "")),
                direction=str(row.get("Direction", "")),
                target_mean_meth=float(row.get("Target_Mean_Meth", 0) or 0),
                background_mean_meth=float(row.get("Background_Mean_Meth", 0) or 0),
                delta_means=float(row.get("Delta_Means", 0) or 0),
                delta_quants=float(row.get("Delta_Quants", 0) or 0),
                delta_maxmin=float(row.get("Delta_MaxMin", 0) or 0),
                ttest_pval=str(row.get("Ttest_Pval", "")),
                mwtest_pval=str(row.get("MWtest_Pval", "")),
                mvalue_ttest=str(row.get("Mvalue_Ttest", "")),
                cleanliness_score=float(row.get("Cleanliness_Score", 0) or 0),
                target_celltype=ct_id,
                cpg_sites=cpg_sites,
            )
            blocks.append(block)

    return blocks


# ---------------------------------------------------------------------------
# Block_Summary format loader (original v2 design)
# ---------------------------------------------------------------------------
def _load_block_summary(xlsx_path: str, cell_type_id: Optional[str] = None) -> List[DMRBlock]:
    """Load DMR blocks from a Block_Summary-format Excel file."""
    block_df = pd.read_excel(xlsx_path, sheet_name="Block_Summary")
    if cell_type_id:
        block_df = block_df[block_df["cell_type_id"] == cell_type_id]

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
            delta_quants=float(row.get("Delta_Quants", 0) or 0),
            delta_maxmin=float(row.get("Delta_MaxMin", 0) or 0),
            ttest_pval=str(row.get("Ttest_Pval", "")),
            mwtest_pval=str(row.get("MWtest_Pval", "")),
            mvalue_ttest=str(row.get("Mvalue_Ttest", "")),
            cleanliness_score=float(row.get("cleanliness_score", 0) or 0),
            target_celltype=str(row.get("Target_CellType", "")),
        )
        blocks.append(block)

    # Load per-CpG data from PerCpG sheets
    sheet_name = f"PerCpG_{cell_type_id}" if cell_type_id else None

    if cell_type_id:
        try:
            percpg_df = pd.read_excel(xlsx_path, sheet_name=f"PerCpG_{cell_type_id}")
        except Exception:
            percpg_df = None
    else:
        all_percpg = []
        xl = pd.ExcelFile(xlsx_path)
        for sheet in xl.sheet_names:
            if sheet.startswith("PerCpG_"):
                df = pd.read_excel(xlsx_path, sheet_name=sheet)
                all_percpg.append(df)
        percpg_df = pd.concat(all_percpg, ignore_index=True) if all_percpg else None

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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def load_dmr_blocks(xlsx_path: str, cell_type_id: Optional[str] = None) -> List[DMRBlock]:
    """
    Load DMR blocks from an Excel file.

    Auto-detects the format (v7.9 or Block_Summary).

    Args:
        xlsx_path: Path to the DMR Excel file
        cell_type_id: If specified, only load blocks for this cell type.
                      Accepts both short IDs (MONO, BCELL) and sheet names (Mono, B).

    Returns:
        List of DMRBlock objects with per-CpG data
    """
    fmt = _detect_format(xlsx_path)

    # Normalize cell_type_id: accept both short IDs and sheet names
    if cell_type_id and cell_type_id in SHEET_TO_CTID:
        # User passed a sheet name (e.g. "Mono") — convert to short ID
        cell_type_id = SHEET_TO_CTID[cell_type_id]

    if fmt == "v79":
        return _load_v79(xlsx_path, cell_type_id)
    elif fmt == "block_summary":
        return _load_block_summary(xlsx_path, cell_type_id)
    else:
        raise ValueError(
            f"Could not detect Excel format in {xlsx_path}. "
            f"Expected either v7.9 format (sheets per cell type) or "
            f"Block_Summary format."
        )


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
