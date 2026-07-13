#!/usr/bin/env python3
"""
Export DMR regions with per-CpG methylation to Excel.

Reads dmr_blocks.json files from pipeline output directories and produces
a multi-sheet Excel workbook with:
  - One sheet per cell type: top N DMR blocks with per-CpG detail
  - A summary sheet: all blocks across all cell types

Usage:
    # Export from a single cell type:
    python export_dmr_excel.py results/MONO/dmr_blocks.json --top 200

    # Export from all 7 cell types (one sheet each):
    python export_dmr_excel.py results/*/dmr_blocks.json --top 200

    # Custom output path:
    python export_dmr_excel.py results/*/dmr_blocks.json --top 200 -o dmr_regions.xlsx
"""

import sys
import json
import glob
import argparse
from pathlib import Path

try:
    import openpyxl
    from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
    from openpyxl.utils import get_column_letter
except ImportError:
    print("Error: openpyxl is required. Install with: pip install openpyxl")
    sys.exit(1)


# ---------------------------------------------------------------------------
# Styling
# ---------------------------------------------------------------------------
HEADER_FONT = Font(bold=True, size=11, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="4472C4", end_color="4472C4", fill_type="solid")
HEADER_ALIGN = Alignment(horizontal="center", vertical="center", wrap_text=True)
THIN_BORDER = Border(
    left=Side(style="thin"), right=Side(style="thin"),
    top=Side(style="thin"), bottom=Side(style="thin"),
)

# Block-level columns (one row per DMR block)
BLOCK_COLUMNS = [
    ("Cell type", "cell_type_id"),
    ("Rank", "rank"),
    ("Seq ID", "seq_id"),
    ("Chromosome", "chrom"),
    ("Start", "start"),
    ("End", "end"),
    ("Block length (bp)", "block_len"),
    ("# CpGs", "num_cpgs"),
    ("Gene", "gene"),
    ("Annotation", "annotation"),
    ("Target mean methylation", "target_mean_meth"),
    ("Background mean methylation", "background_mean_meth"),
    ("Delta means", "delta_means"),
    ("Cleanliness score", "cleanliness_score"),
]

# Per-CpG columns (one row per CpG within each block)
CPG_COLUMNS = [
    ("Cell type", "cell_type_id"),
    ("Seq ID", "seq_id"),
    ("Rank", "rank"),
    ("Chromosome", "chrom"),
    ("Block start", "start"),
    ("Block end", "end"),
    ("Gene", "gene"),
    ("CpG label", "label"),
    ("CpG position", "position"),
    ("CpG index", "global_idx"),
    ("Target mean beta", "target_mean_beta"),
    ("Background mean beta", "background_mean_beta"),
    ("Delta beta", "delta_beta"),
]


def load_blocks(json_path: str) -> list:
    """Load DMR blocks from a JSON file."""
    with open(json_path, "r") as f:
        return json.load(f)


def format_block_row(block: dict) -> dict:
    """Extract block-level columns from a block dict."""
    row = {}
    for display_name, key in BLOCK_COLUMNS:
        row[display_name] = block.get(key, "")
    return row


def format_cpg_rows(block: dict) -> list:
    """Extract one row per CpG from a block dict."""
    rows = []
    for cpg in block.get("cpg_sites", []):
        row = {}
        for display_name, key in CPG_COLUMNS:
            if key in cpg:
                row[display_name] = cpg[key]
            else:
                row[display_name] = block.get(key, "")
        # Compute delta_beta if not present
        if "Delta beta" not in row or row["Delta beta"] == "":
            tg = cpg.get("target_mean_beta", 0)
            bg = cpg.get("background_mean_beta", 0)
            row["Delta beta"] = round(abs(tg - bg), 4)
        rows.append(row)
    return rows


def write_sheet(ws, headers, rows, freeze="A2"):
    """Write headers and rows to a worksheet with styling."""
    # Headers
    for col_idx, header in enumerate(headers, 1):
        cell = ws.cell(row=1, column=col_idx, value=header)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = HEADER_ALIGN
        cell.border = THIN_BORDER

    # Data rows
    for row_idx, row_data in enumerate(rows, 2):
        for col_idx, header in enumerate(headers, 1):
            cell = ws.cell(row=row_idx, column=col_idx, value=row_data.get(header, ""))
            cell.border = THIN_BORDER

    # Auto-width (approximate)
    for col_idx, header in enumerate(headers, 1):
        max_len = len(str(header))
        for row_idx in range(2, min(len(rows) + 2, 200)):  # sample first 200 rows
            val = ws.cell(row=row_idx, column=col_idx).value
            if val is not None:
                max_len = max(max_len, len(str(val)))
        ws.column_dimensions[get_column_letter(col_idx)].width = min(max_len + 3, 35)

    # Freeze header row
    ws.freeze_panes = freeze


def export_excel(json_paths, top_n, output_path):
    """Export DMR blocks from multiple JSON files to a multi-sheet Excel."""
    wb = openpyxl.Workbook()

    # Remove default sheet
    wb.remove(wb.active)

    all_block_rows = []
    all_cpg_rows = []

    for json_path in sorted(json_paths):
        blocks = load_blocks(json_path)
        if not blocks:
            print(f"  Warning: no blocks in {json_path}, skipping")
            continue

        # Sort by rank (rank 1 = best), take top N
        blocks_sorted = sorted(blocks, key=lambda b: b.get("rank", 999))
        top_blocks = blocks_sorted[:top_n]

        cell_type = top_blocks[0].get("cell_type_id", "Unknown")
        print(f"  {json_path}: {len(blocks)} blocks → top {len(top_blocks)} for {cell_type}")

        # Per-cell-type sheet: block-level summary
        block_rows = [format_block_row(b) for b in top_blocks]
        block_headers = [name for name, _ in BLOCK_COLUMNS]
        ws = wb.create_sheet(title=cell_type[:31])  # Excel sheet name max 31 chars
        write_sheet(ws, block_headers, block_rows)

        # Collect for summary sheets
        all_block_rows.extend(block_rows)
        for b in top_blocks:
            all_cpg_rows.extend(format_cpg_rows(b))

    # Summary sheet: all blocks
    if all_block_rows:
        block_headers = [name for name, _ in BLOCK_COLUMNS]
        ws_summary = wb.create_sheet(title="All blocks summary", index=0)
        write_sheet(ws_summary, block_headers, all_block_rows)

    # Per-CpG sheet: all CpGs across all cell types
    if all_cpg_rows:
        cpg_headers = [name for name, _ in CPG_COLUMNS]
        ws_cpg = wb.create_sheet(title="Per-CpG detail")
        write_sheet(ws_cpg, cpg_headers, all_cpg_rows)

    wb.save(output_path)
    print(f"\nSaved: {output_path}")
    print(f"  {len(all_block_rows)} blocks across {len(json_paths)} cell type(s)")
    print(f"  {len(all_cpg_rows)} CpG sites total")


def main():
    parser = argparse.ArgumentParser(
        description="Export DMR regions with per-CpG methylation to Excel."
    )
    parser.add_argument(
        "inputs", nargs="+",
        help="Path(s) to dmr_blocks.json file(s). Supports glob patterns."
    )
    parser.add_argument(
        "--top", type=int, default=200,
        help="Number of top DMR regions per cell type (default: 200)"
    )
    parser.add_argument(
        "-o", "--output", default="dmr_regions.xlsx",
        help="Output Excel file path (default: dmr_regions.xlsx)"
    )
    args = parser.parse_args()

    # Expand glob patterns
    json_paths = []
    for pattern in args.inputs:
        matches = sorted(glob.glob(pattern))
        if matches:
            json_paths.extend(matches)
        elif Path(pattern).exists():
            json_paths.append(pattern)
        else:
            print(f"Warning: no files matched '{pattern}'")

    if not json_paths:
        print("Error: no dmr_blocks.json files found.")
        sys.exit(1)

    # Deduplicate
    json_paths = sorted(set(json_paths))

    print(f"Exporting top {args.top} DMR regions from {len(json_paths)} file(s):")
    export_excel(json_paths, args.top, args.output)


if __name__ == "__main__":
    main()
