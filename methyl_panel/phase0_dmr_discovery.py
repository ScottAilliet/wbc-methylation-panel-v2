#!/usr/bin/env python3
"""
Phase 0: DMR Discovery via wgbstools find_markers
=================================================

Discovers differentially methylated regions (DMRs) from raw beta files
using wgbstools find_markers, then extracts per-CpG methylation data
and computes cleanliness scores.

This module makes the pipeline fully self-contained: no pre-computed
DMR Excel file is needed. Anyone with the downloaded beta files and
blocks file can discover DMRs from scratch.

Workflow:
    1. Generate a groups CSV (target vs background) for the requested cell type
    2. Generate a beta list file (paths to all .beta files)
    3. Run wgbstools find_markers via subprocess
    4. Parse the output BED file into DMRBlock objects
    5. Extract per-CpG methylation from beta files
    6. Compute cleanliness scores
    7. Save as JSON (same format as step 1 output)

Usage (called by pipeline.py when --discover-dmrs is set):
    from methyl_panel.phase0_dmr_discovery import discover_dmrs
    blocks = discover_dmrs(cell_type_id, beta_dir, blocks_file, groups_csv, out_dir)
"""

import os
import sys
import glob
import json
import shutil
import subprocess
import warnings
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Tuple

# Suppress harmless RuntimeWarnings from numpy (empty slices when all NaN)
warnings.filterwarnings("ignore", category=RuntimeWarning)

from methyl_panel.phase1_dmr_loader import DMRBlock, CpGSite


# ---------------------------------------------------------------------------
# Cell type → target groups mapping
# ---------------------------------------------------------------------------
# Each cell type ID maps to a list of atlas group names that constitute
# the target population. Background is restricted to other blood cell
# types only — this is a WBC assay, so we need to distinguish monocytes
# from other blood cells, not from liver or brain tissue.
CELL_TYPE_TARGETS: Dict[str, List[str]] = {
    "MONO":  ["Blood-Monocytes"],
    "BCELL": ["Blood-B", "Blood-B-Mem"],
    "NK":    ["Blood-NK"],
    "GRAN":  ["Blood-Granulocytes"],
    "CD3T":  [
        "Blood-T-CD3", "Blood-T-CD4", "Blood-T-CD8",
        "Blood-T-CenMem-CD4", "Blood-T-EffMem-CD4", "Blood-T-EffMem-CD8",
        "Blood-T-Eff-CD8", "Blood-T-Naive-CD4", "Blood-T-Naive-CD8",
    ],
    "CD4T":  [
        "Blood-T-CD4", "Blood-T-CenMem-CD4",
        "Blood-T-EffMem-CD4", "Blood-T-Naive-CD4",
    ],
    "CD8T":  [
        "Blood-T-CD8", "Blood-T-Eff-CD8",
        "Blood-T-EffMem-CD8", "Blood-T-Naive-CD8",
    ],
}

# All blood cell type groups in the atlas (14 groups, 36 samples total).
# Background for any WBC assay is restricted to these groups minus the
# target groups for the cell type being assayed.
BLOOD_GROUPS: List[str] = [
    "Blood-B", "Blood-B-Mem", "Blood-Granulocytes", "Blood-Monocytes",
    "Blood-NK", "Blood-T-CD3", "Blood-T-CD4", "Blood-T-CD8",
    "Blood-T-CenMem-CD4", "Blood-T-Eff-CD8", "Blood-T-EffMem-CD4",
    "Blood-T-EffMem-CD8", "Blood-T-Naive-CD4", "Blood-T-Naive-CD8",
]

# All valid cell type IDs
VALID_CELL_TYPES = list(CELL_TYPE_TARGETS.keys())


# ---------------------------------------------------------------------------
# Groups file generation
# ---------------------------------------------------------------------------
def generate_groups_file(
    cell_type_id: str,
    atlas_groups_csv: str,
    output_path: str,
    beta_dir: Optional[str] = None,
) -> str:
    """
    Generate a wgbstools-compatible groups CSV for one cell type.

    Target samples (all subtypes merged) get group = cell_type_id.
    Background is restricted to other blood cell types only — this is
    a WBC assay, so we distinguish the target from other blood cells,
    not from non-blood tissues.

    If beta_dir is provided, only samples that have a corresponding
    .beta file in that directory are included in the groups file.
    This prevents find_markers from failing on missing beta files.

    Args:
        cell_type_id: One of MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T
        atlas_groups_csv: Path to data/full_atlas_groups.csv
        output_path: Where to write the groups CSV
        beta_dir: If provided, filter to samples with .beta files here

    Returns:
        Path to the generated groups CSV
    """
    if cell_type_id not in CELL_TYPE_TARGETS:
        raise ValueError(
            f"Unknown cell type '{cell_type_id}'. "
            f"Valid options: {VALID_CELL_TYPES}"
        )

    target_groups = CELL_TYPE_TARGETS[cell_type_id]
    df = pd.read_csv(atlas_groups_csv)

    # Restrict to blood cell types only (exclude non-blood tissues)
    df = df[df['group'].isin(BLOOD_GROUPS)].copy()

    # The atlas groups CSV has columns: name, group
    # wgbstools renames the first column to 'fname', so 'name' is fine
    # Merge target subtypes into one group
    df['group'] = df['group'].apply(
        lambda g: cell_type_id if g in target_groups else 'background'
    )

    # Filter to samples that have beta files, if beta_dir is provided
    if beta_dir:
        import glob
        available = set()
        for bp in glob.glob(os.path.join(beta_dir, "*.beta")):
            stem = os.path.basename(bp)
            if stem.endswith('.beta'):
                stem = stem[:-5]
            available.add(stem)
        name_col = df.columns[0]
        before = len(df)
        df = df[df[name_col].isin(available)]
        after = len(df)
        n_target = (df['group'] != 'background').sum()
        n_bg = (df['group'] == 'background').sum()
        print(f"    Groups file: {after}/{before} samples have beta files "
              f"({n_target} target, {n_bg} background)")

    # Write with 'fname' column name for clarity (wgbstools renames first col anyway)
    df.to_csv(output_path, index=False)
    return output_path


# ---------------------------------------------------------------------------
# Beta list file generation
# ---------------------------------------------------------------------------
def generate_beta_list(beta_dir: str, output_path: str) -> str:
    """
    Generate a text file listing all .beta file paths.

    Args:
        beta_dir: Directory containing .beta files
        output_path: Where to write the beta list

    Returns:
        Path to the generated beta list
    """
    beta_files = sorted(glob.glob(os.path.join(beta_dir, "*.beta")))
    if not beta_files:
        raise FileNotFoundError(
            f"No .beta files found in {beta_dir}. "
            f"Run download_data.sh first."
        )
    with open(output_path, 'w') as f:
        for path in beta_files:
            f.write(path + '\n')
    return output_path


# ---------------------------------------------------------------------------
# Run wgbstools find_markers
# ---------------------------------------------------------------------------
def find_wgbstools(wgbstools_path: Optional[str] = None) -> str:
    """
    Locate the wgbstools executable.

    Priority:
        1. Explicit path provided
        2. 'wgbstools' on PATH
        3. ./wgbs_tools/wgbstools (cloned by install_macos.sh)
    """
    if wgbstools_path and os.path.isfile(wgbstools_path):
        return wgbstools_path

    on_path = shutil.which('wgbstools')
    if on_path:
        return on_path

    local = os.path.join(os.getcwd(), 'wgbs_tools', 'wgbstools')
    if os.path.isfile(local):
        return local

    raise FileNotFoundError(
        "wgbstools not found. Install it with: "
        "git clone https://github.com/nloyfer/wgbs_tools.git && "
        "cd wgbs_tools && python3 setup.py && cd .. && "
        "ln -sf \"$(pwd)/wgbs_tools/wgbstools\" .venv/bin/wgbstools"
    )


def run_find_markers(
    groups_file: str,
    beta_list: str,
    blocks_file: str,
    out_dir: str,
    cell_type_id: str,
    threads: int = 2,
    top_n: int = 200,
    delta_means: float = 0.3,
    min_cpg: int = 3,
    wgbstools_path: Optional[str] = None,
) -> str:
    """
    Run wgbstools find_markers for one cell type.

    Args:
        groups_file: Path to groups CSV
        beta_list: Path to beta list text file
        blocks_file: Path to blocks BED file
        out_dir: Output directory for find_markers
        cell_type_id: Target group name in the groups file
        threads: Number of threads
        top_n: Number of top markers to keep
        delta_means: Minimum methylation difference
        min_cpg: Minimum CpGs per block
        wgbstools_path: Path to wgbstools executable

    Returns:
        Path to the output BED file (Markers.{cell_type_id}.bed)
    """
    wgbstools = find_wgbstools(wgbstools_path)
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        wgbstools, 'find_markers',
        '--blocks_path', blocks_file,
        '--groups_file', groups_file,
        '--beta_list_file', beta_list,
        '--targets', cell_type_id,
        '--only_hypo',
        '--delta_means', str(delta_means),
        '--min_cpg', str(min_cpg),
        '--sort_by', 'delta_means',
        '--top', str(top_n),
        '--header',
        '--threads', str(threads),
        '--out_dir', out_dir,
    ]

    print(f"  Running: {' '.join(cmd[:6])} ... --out_dir {out_dir}")
    result = subprocess.run(
        cmd, capture_output=True, text=True, timeout=3600
    )

    if result.returncode != 0:
        # Print stderr for debugging
        if result.stderr:
            print(f"  wgbstools stderr: {result.stderr[:500]}")
        raise RuntimeError(
            f"wgbstools find_markers failed with return code {result.returncode}"
        )

    if result.stderr:
        # find_markers prints progress to stderr
        for line in result.stderr.strip().split('\n')[-5:]:
            print(f"  {line}")

    bed_path = os.path.join(out_dir, f'Markers.{cell_type_id}.bed')
    if not os.path.isfile(bed_path):
        raise FileNotFoundError(
            f"Expected output not found: {bed_path}"
        )

    return bed_path


# ---------------------------------------------------------------------------
# Parse find_markers BED output
# ---------------------------------------------------------------------------
# Column names from find_markers dump_results:
MARKERS_BED_COLUMNS = [
    'chr', 'start', 'end', 'startCpG', 'endCpG',
    'target', 'region', 'lenCpG', 'bp',
    'tg_mean', 'bg_mean', 'delta_means', 'delta_quants',
    'delta_maxmin', 'ttest', 'mw_test', 'mvalue_ttest',
    'direction', 'anno', 'gene',
]


def parse_markers_bed(bed_path: str, cell_type_id: str) -> List[DMRBlock]:
    """
    Parse a find_markers BED file into DMRBlock objects.

    The BED file has comment lines starting with '#' (sample lists)
    followed by a header line and data lines.

    Args:
        bed_path: Path to Markers.{cell_type_id}.bed
        cell_type_id: Cell type ID for these markers

    Returns:
        List of DMRBlock objects (without per-CpG data yet)
    """
    # Read the file, skipping comment lines (but keeping the #chr header)
    with open(bed_path) as f:
        lines = [l for l in f if not l.startswith('#>') and not l.startswith('#<')]

    if not lines:
        return []

    # Parse with pandas
    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    # Rename #chr to chr if present
    if '#chr' in df.columns:
        df = df.rename(columns={'#chr': 'chr'})

    blocks = []
    for rank, (_, row) in enumerate(df.iterrows(), 1):
        chrom = str(row['chr'])
        start = int(row['start'])
        end = int(row['end'])
        start_cpg = int(row['startCpG'])
        end_cpg = int(row['endCpG'])
        num_cpgs = end_cpg - start_cpg

        # Parse lenCpG string like "7CpGs"
        len_cpg_str = str(row.get('lenCpG', ''))
        if len_cpg_str.endswith('CpGs'):
            try:
                num_cpgs = int(len_cpg_str.replace('CpGs', ''))
            except ValueError:
                pass

        # Parse bp string like "120bp"
        bp_str = str(row.get('bp', ''))
        block_len = end - start
        if bp_str.endswith('bp'):
            try:
                block_len = int(bp_str.replace('bp', ''))
            except ValueError:
                pass

        seq_id = f"{cell_type_id}_{rank:04d}"

        block = DMRBlock(
            cell_type_id=cell_type_id,
            seq_id=seq_id,
            rank=rank,
            chrom=chrom,
            start=start,
            end=end,
            num_cpgs=num_cpgs,
            block_len=block_len,
            gene=str(row.get('gene', '')) if pd.notna(row.get('gene', '')) else '',
            annotation=str(row.get('anno', '')) if pd.notna(row.get('anno', '')) else '',
            direction=str(row.get('direction', '')) if pd.notna(row.get('direction', '')) else '',
            target_mean_meth=float(row.get('tg_mean', 0) or 0),
            background_mean_meth=float(row.get('bg_mean', 0) or 0),
            delta_means=float(row.get('delta_means', 0) or 0),
            delta_quants=float(row.get('delta_quants', 0) or 0),
            delta_maxmin=float(row.get('delta_maxmin', 0) or 0),
            ttest_pval=str(row.get('ttest', '')) if pd.notna(row.get('ttest', '')) else '',
            mwtest_pval=str(row.get('mw_test', '')) if pd.notna(row.get('mw_test', '')) else '',
            mvalue_ttest=str(row.get('mvalue_ttest', '')) if pd.notna(row.get('mvalue_ttest', '')) else '',
            cleanliness_score=0.0,  # Will be computed from per-CpG data
            target_celltype=cell_type_id,
            cpg_sites=[],  # Will be filled from beta files
        )
        blocks.append(block)

    return blocks


# ---------------------------------------------------------------------------
# Per-CpG methylation extraction from beta files
# ---------------------------------------------------------------------------
def read_beta_region(
    beta_path: str,
    start_cpg: int,
    end_cpg: int,
    min_cov: int = 5,
) -> np.ndarray:
    """
    Read per-CpG methylation fractions from a wgbstools beta file.

    Beta file format:
        - Binary, 2 x uint8 per CpG: (#meth, #covered)
        - CpG index is 1-based (global)
        - Byte offset = (start_cpg - 1) * 2

    Args:
        beta_path: Path to .beta file
        start_cpg: Start CpG index (1-based, inclusive)
        end_cpg: End CpG index (1-based, exclusive)
        min_cov: Minimum coverage to include

    Returns:
        np.array of float, shape (end_cpg - start_cpg,), NaN where cov < min_cov
    """
    n_cpgs = end_cpg - start_cpg
    if n_cpgs < 1:
        return np.array([])

    offset = (start_cpg - 1) * 2
    try:
        with open(beta_path, "rb") as f:
            f.seek(offset)
            raw = np.frombuffer(
                f.read(n_cpgs * 2), dtype=np.uint8
            ).reshape(-1, 2)
        meth = raw[:, 0].astype(float)
        cov = raw[:, 1].astype(float)
        fracs = np.where(cov >= min_cov, meth / np.maximum(cov, 1), np.nan)
        return fracs
    except Exception as e:
        print(f"    WARNING: could not read {os.path.basename(beta_path)}: {e}")
        return np.full(n_cpgs, np.nan)


def _match_beta_files(
    groups_csv: str,
    beta_dir: str,
) -> Tuple[List[str], List[str]]:
    """
    Match target and background sample names to beta file paths.

    Args:
        groups_csv: Path to the generated groups CSV (target/background)
        beta_dir: Directory containing .beta files

    Returns:
        (target_betas, background_betas) — lists of file paths
    """
    df = pd.read_csv(groups_csv)
    # First column is sample name (renamed to 'fname' by wgbstools, but
    # we wrote it with original column name)
    name_col = df.columns[0]

    all_betas = glob.glob(os.path.join(beta_dir, "*.beta"))
    # Build a lookup: sample_name -> beta_path
    beta_lookup = {}
    for bp in all_betas:
        basename = os.path.basename(bp)
        # Remove .beta extension
        stem = basename[:-5] if basename.endswith('.beta') else basename
        beta_lookup[stem] = bp

    target_names = df[df['group'] != 'background'][name_col].tolist()
    bg_names = df[df['group'] == 'background'][name_col].tolist()

    target_betas = []
    for name in target_names:
        if name in beta_lookup:
            target_betas.append(beta_lookup[name])
        else:
            print(f"    WARNING: no beta file for target sample '{name}'")

    bg_betas = []
    for name in bg_names:
        if name in beta_lookup:
            bg_betas.append(beta_lookup[name])
        # Don't warn for missing background — there may be many

    return target_betas, bg_betas


def extract_percpg_methylation(
    blocks: List[DMRBlock],
    target_betas: List[str],
    background_betas: List[str],
    min_cov: int = 5,
    max_bg_samples: int = 30,
) -> None:
    """
    Extract per-CpG methylation for each block and populate cpg_sites.

    For each block, reads target and background beta files at the block's
    CpG indices (startCpG to endCpG-1), computes per-CpG mean methylation,
    and creates CpGSite objects.

    Also computes the cleanliness score.

    Args:
        blocks: List of DMRBlock objects (modified in place)
        target_betas: List of target beta file paths
        background_betas: List of background beta file paths
        min_cov: Minimum coverage per CpG site
        max_bg_samples: Maximum background samples to use (for speed)
    """
    # Subsample background if too many
    bg_betas = background_betas[:max_bg_samples]
    print(f"  Per-CpG extraction: {len(target_betas)} target, "
          f"{len(bg_betas)} background samples")

    for i, block in enumerate(blocks):
        start_cpg = None
        end_cpg = None

        # We need the CpG indices. These were in the BED file but not
        # stored in DMRBlock. We'll re-derive them from the BED file
        # or use the block's start/end genomic positions.
        # Actually, we stored them during parsing — let's use a workaround:
        # The DMRBlock doesn't have startCpG/endCpG fields, so we need to
        # pass them separately. We'll use the block's num_cpgs and
        # distribute CpGs evenly if we can't get indices.

        # For now, skip per-CpG extraction if we don't have CpG indices
        # The cleanliness score will be 0, and CpG positions will be
        # found by sequence scanning in step 2.

        # Create placeholder CpG sites with even distribution
        if block.num_cpgs > 0 and not block.cpg_sites:
            step = max(1, (block.end - block.start) // max(1, block.num_cpgs))
            for j in range(block.num_cpgs):
                pos = block.start + j * step
                block.cpg_sites.append(CpGSite(
                    label=f"CpG_{j+1}",
                    position=pos,
                    global_idx=pos,  # Placeholder; will be corrected in step 2
                    target_mean_beta=block.target_mean_meth,
                    background_mean_beta=block.background_mean_meth,
                    delta_beta=abs(block.target_mean_meth - block.background_mean_meth),
                ))

        # Compute cleanliness score from block-level stats
        # (simplified: use delta_means as a proxy)
        block.cleanliness_score = min(1.0, block.delta_means / 0.5)

        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(blocks)} blocks")


def extract_percpg_methylation_with_indices(
    blocks: List[DMRBlock],
    cpg_indices: List[Tuple[int, int]],  # (startCpG, endCpG) per block
    target_betas: List[str],
    background_betas: List[str],
    min_cov: int = 5,
    max_bg_samples: int = 30,
) -> None:
    """
    Extract per-CpG methylation using CpG indices from find_markers.

    This is the full version that reads beta files at the correct CpG
    indices and computes proper cleanliness scores.

    Args:
        blocks: List of DMRBlock objects (modified in place)
        cpg_indices: List of (startCpG, endCpG) tuples, one per block
        target_betas: List of target beta file paths
        background_betas: List of background beta file paths
        min_cov: Minimum coverage per CpG site
        max_bg_samples: Maximum background samples to use (for speed)
    """
    bg_betas = background_betas[:max_bg_samples]
    print(f"  Per-CpG extraction: {len(target_betas)} target, "
          f"{len(bg_betas)} background samples")

    for i, (block, (start_cpg, end_cpg)) in enumerate(zip(blocks, cpg_indices)):
        n_cpgs = end_cpg - start_cpg
        if n_cpgs < 1:
            continue

        # Read per-CpG methylation from all target and background beta files
        tg_mat = np.array([
            read_beta_region(bp, start_cpg, end_cpg, min_cov)
            for bp in target_betas
        ]) if target_betas else np.full((1, n_cpgs), np.nan)

        bg_mat = np.array([
            read_beta_region(bp, start_cpg, end_cpg, min_cov)
            for bp in bg_betas
        ]) if bg_betas else np.full((1, n_cpgs), np.nan)

        # Per-CpG mean methylation
        tg_cpg_mean = np.nanmean(tg_mat, axis=0)
        bg_cpg_mean = np.nanmean(bg_mat, axis=0)

        # Compute cleanliness score
        score = _cleanliness_score(tg_mat, bg_mat)

        # Build CpG sites
        # Genomic positions: distribute evenly across the block
        # (precise positions will be found by sequence scanning in step 2)
        block.cpg_sites = []
        step = max(1, (block.end - block.start) // max(1, n_cpgs))
        for j in range(n_cpgs):
            pos = block.start + j * step
            tg_val = float(tg_cpg_mean[j]) if not np.isnan(tg_cpg_mean[j]) else 0.0
            bg_val = float(bg_cpg_mean[j]) if not np.isnan(bg_cpg_mean[j]) else 0.0
            block.cpg_sites.append(CpGSite(
                label=f"CpG_{j+1}",
                position=pos,
                global_idx=start_cpg + j,  # Actual CpG index
                target_mean_beta=tg_val,
                background_mean_beta=bg_val,
                delta_beta=abs(tg_val - bg_val),
            ))

        block.cleanliness_score = score

        if (i + 1) % 50 == 0:
            print(f"  Processed {i+1}/{len(blocks)} blocks")


# ---------------------------------------------------------------------------
# Cleanliness score (from repo 1's per_cpg_marker_analysis.py)
# ---------------------------------------------------------------------------
def _cleanliness_score(
    tg_matrix: np.ndarray,
    bg_matrix: np.ndarray,
) -> float:
    """
    Score how suitable a marker block is for bisulfite PCR.

    Components (each 0-1, higher = better):
        A. tg_near_zero: 1 - max(per-CpG mean in target)
        B. bg_near_one: min(per-CpG mean in background)
        C. tg_consistency: 1 - 2*mean(per-CpG std across target samples)
        D. bg_consistency: 1 - 2*mean(per-CpG std across background samples)
        E. coverage: fraction of CpGs with data in >=50% of samples

    Final score = mean(A, B, C, D) * E

    Args:
        tg_matrix: (n_target_samples, n_cpgs) methylation fractions
        bg_matrix: (n_bg_samples, n_cpgs) methylation fractions

    Returns:
        Cleanliness score (0-1)
    """
    tg_cpg_mean = np.nanmean(tg_matrix, axis=0)
    bg_cpg_mean = np.nanmean(bg_matrix, axis=0)

    A = max(0.0, float(1.0 - np.nanmax(tg_cpg_mean))) if not np.all(np.isnan(tg_cpg_mean)) else 0.0
    B = max(0.0, float(np.nanmin(bg_cpg_mean))) if not np.all(np.isnan(bg_cpg_mean)) else 0.0
    C = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(tg_matrix, axis=0))))
    D = max(0.0, 1.0 - 2.0 * float(np.nanmean(np.nanstd(bg_matrix, axis=0))))

    tg_cov = float(np.mean(np.mean(~np.isnan(tg_matrix), axis=0) >= 0.5))
    bg_cov = float(np.mean(np.mean(~np.isnan(bg_matrix), axis=0) >= 0.5))
    E = (tg_cov + bg_cov) / 2.0

    score = np.mean([A, B, C, D]) * E
    return round(float(score), 4)


# ---------------------------------------------------------------------------
# Parse BED with CpG indices
# ---------------------------------------------------------------------------
def parse_markers_bed_with_indices(
    bed_path: str,
    cell_type_id: str,
) -> Tuple[List[DMRBlock], List[Tuple[int, int]]]:
    """
    Parse a find_markers BED file, returning both DMRBlocks and CpG indices.

    Args:
        bed_path: Path to Markers.{cell_type_id}.bed
        cell_type_id: Cell type ID for these markers

    Returns:
        (blocks, cpg_indices) where cpg_indices[i] = (startCpG, endCpG) for block i
    """
    with open(bed_path) as f:
        lines = [l for l in f if not l.startswith('#>') and not l.startswith('#<')]

    if not lines:
        return [], []

    from io import StringIO
    df = pd.read_csv(StringIO(''.join(lines)), sep='\t')

    if '#chr' in df.columns:
        df = df.rename(columns={'#chr': 'chr'})

    blocks = []
    cpg_indices = []

    for rank, (_, row) in enumerate(df.iterrows(), 1):
        chrom = str(row['chr'])
        start = int(row['start'])
        end = int(row['end'])
        start_cpg = int(row['startCpG'])
        end_cpg = int(row['endCpG'])
        num_cpgs = end_cpg - start_cpg

        len_cpg_str = str(row.get('lenCpG', ''))
        if len_cpg_str.endswith('CpGs'):
            try:
                num_cpgs = int(len_cpg_str.replace('CpGs', ''))
            except ValueError:
                pass

        bp_str = str(row.get('bp', ''))
        block_len = end - start
        if bp_str.endswith('bp'):
            try:
                block_len = int(bp_str.replace('bp', ''))
            except ValueError:
                pass

        seq_id = f"{cell_type_id}_{rank:04d}"

        block = DMRBlock(
            cell_type_id=cell_type_id,
            seq_id=seq_id,
            rank=rank,
            chrom=chrom,
            start=start,
            end=end,
            num_cpgs=num_cpgs,
            block_len=block_len,
            gene=str(row.get('gene', '')) if pd.notna(row.get('gene', '')) else '',
            annotation=str(row.get('anno', '')) if pd.notna(row.get('anno', '')) else '',
            direction=str(row.get('direction', '')) if pd.notna(row.get('direction', '')) else '',
            target_mean_meth=float(row.get('tg_mean', 0) or 0),
            background_mean_meth=float(row.get('bg_mean', 0) or 0),
            delta_means=float(row.get('delta_means', 0) or 0),
            delta_quants=float(row.get('delta_quants', 0) or 0),
            delta_maxmin=float(row.get('delta_maxmin', 0) or 0),
            ttest_pval=str(row.get('ttest', '')) if pd.notna(row.get('ttest', '')) else '',
            mwtest_pval=str(row.get('mw_test', '')) if pd.notna(row.get('mw_test', '')) else '',
            mvalue_ttest=str(row.get('mvalue_ttest', '')) if pd.notna(row.get('mvalue_ttest', '')) else '',
            cleanliness_score=0.0,
            target_celltype=cell_type_id,
            cpg_sites=[],
        )
        blocks.append(block)
        cpg_indices.append((start_cpg, end_cpg))

    return blocks, cpg_indices


# ---------------------------------------------------------------------------
# Save blocks as JSON (same format as step 1)
# ---------------------------------------------------------------------------
def save_blocks_json(blocks: List[DMRBlock], output_path: str) -> None:
    """Save DMRBlocks as JSON in the same format as step 1 output."""
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

    with open(output_path, 'w') as f:
        json.dump(blocks_data, f, indent=2)
    print(f"  Saved {len(blocks)} blocks to {output_path}")


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------
def discover_dmrs(
    cell_type_id: str,
    beta_dir: str,
    blocks_file: str,
    groups_csv: str,
    out_dir: str,
    threads: int = 2,
    top_n: int = 200,
    delta_means: float = 0.3,
    min_cpg: int = 3,
    wgbstools_path: Optional[str] = None,
    max_bg_samples: int = 30,
    skip_find_markers: bool = False,
) -> List[DMRBlock]:
    """
    Full DMR discovery pipeline for one cell type.

    Args:
        cell_type_id: One of MONO, BCELL, NK, GRAN, CD3T, CD4T, CD8T
        beta_dir: Directory containing .beta files
        blocks_file: Path to blocks BED file
        groups_csv: Path to data/full_atlas_groups.csv
        out_dir: Output directory
        threads: Threads for find_markers
        top_n: Top N markers to keep
        delta_means: Min methylation difference
        min_cpg: Min CpGs per block
        wgbstools_path: Path to wgbstools executable
        max_bg_samples: Max background samples for per-CpG extraction
        skip_find_markers: If True, skip find_markers (use existing BED)

    Returns:
        List of DMRBlock objects with per-CpG data
    """
    os.makedirs(out_dir, exist_ok=True)

    # Step 0a: Generate groups file
    print(f"\n  Step 0a: Generating groups file for {cell_type_id}")
    groups_file = os.path.join(out_dir, f'groups_{cell_type_id}.csv')
    generate_groups_file(cell_type_id, groups_csv, groups_file, beta_dir=beta_dir)
    n_target = sum(1 for _, r in pd.read_csv(groups_file).iterrows()
                   if r['group'] != 'background')
    print(f"    {n_target} target samples, "
          f"{len(pd.read_csv(groups_file)) - n_target} background")

    # Step 0b: Generate beta list
    print(f"  Step 0b: Generating beta list")
    beta_list = os.path.join(out_dir, 'beta_list.txt')
    generate_beta_list(beta_dir, beta_list)
    print(f"    {sum(1 for _ in open(beta_list))} beta files")

    # Step 0c: Run find_markers
    fm_out_dir = os.path.join(out_dir, 'find_markers_output')
    bed_path = os.path.join(fm_out_dir, f'Markers.{cell_type_id}.bed')

    if skip_find_markers and os.path.isfile(bed_path):
        print(f"  Step 0c: Skipping find_markers (using existing {bed_path})")
    else:
        print(f"  Step 0c: Running wgbstools find_markers")
        bed_path = run_find_markers(
            groups_file=groups_file,
            beta_list=beta_list,
            blocks_file=blocks_file,
            out_dir=fm_out_dir,
            cell_type_id=cell_type_id,
            threads=threads,
            top_n=top_n,
            delta_means=delta_means,
            min_cpg=min_cpg,
            wgbstools_path=wgbstools_path,
        )

    # Step 0d: Parse BED output
    print(f"  Step 0d: Parsing find_markers output")
    blocks, cpg_indices = parse_markers_bed_with_indices(bed_path, cell_type_id)
    print(f"    Found {len(blocks)} markers")

    if not blocks:
        print("    WARNING: No markers found! Check find_markers parameters.")
        return []

    # Step 0e: Extract per-CpG methylation
    print(f"  Step 0e: Extracting per-CpG methylation from beta files")
    target_betas, bg_betas = _match_beta_files(groups_file, beta_dir)
    extract_percpg_methylation_with_indices(
        blocks, cpg_indices, target_betas, bg_betas,
        max_bg_samples=max_bg_samples,
    )

    # Report cleanliness scores
    good = sum(1 for b in blocks if b.cleanliness_score >= 0.6)
    print(f"    Cleanliness: {good}/{len(blocks)} blocks score >= 0.6")

    # Step 0f: Save as JSON
    json_path = os.path.join(out_dir, 'dmr_blocks.json')
    save_blocks_json(blocks, json_path)

    return blocks


# ---------------------------------------------------------------------------
# CLI for standalone testing
# ---------------------------------------------------------------------------
def main():
    """CLI entry point for testing phase0 independently."""
    import argparse
    parser = argparse.ArgumentParser(
        description='DMR Discovery via wgbstools find_markers'
    )
    parser.add_argument('--cell-type', required=True,
                        help=f'Cell type ID: {", ".join(VALID_CELL_TYPES)}')
    parser.add_argument('--beta-dir', default='data/beta_files/',
                        help='Directory with .beta files')
    parser.add_argument('--blocks-file',
                        default='data/GSE186458_blocks.s205.bed.gz',
                        help='Path to blocks BED file')
    parser.add_argument('--groups-csv', default='data/full_atlas_groups.csv',
                        help='Path to full atlas groups CSV')
    parser.add_argument('--output-dir', default='dmr_output/',
                        help='Output directory')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--top', type=int, default=200)
    parser.add_argument('--max-bg-samples', type=int, default=30)
    parser.add_argument('--skip-find-markers', action='store_true',
                        help='Skip find_markers (use existing BED)')
    args = parser.parse_args()

    blocks = discover_dmrs(
        cell_type_id=args.cell_type,
        beta_dir=args.beta_dir,
        blocks_file=args.blocks_file,
        groups_csv=args.groups_csv,
        out_dir=args.output_dir,
        threads=args.threads,
        top_n=args.top,
        max_bg_samples=args.max_bg_samples,
        skip_find_markers=args.skip_find_markers,
    )

    print(f"\n=== DMR Discovery complete ===")
    print(f"  Cell type: {args.cell_type}")
    print(f"  Markers found: {len(blocks)}")
    if blocks:
        print(f"\n  Top 5 markers:")
        for b in blocks[:5]:
            print(f"    {b.seq_id}: {b.chrom}:{b.start}-{b.end} "
                  f"{b.num_cpgs}CpGs delta={b.delta_means:.3f} "
                  f"score={b.cleanliness_score:.3f} gene={b.gene}")


if __name__ == "__main__":
    main()
