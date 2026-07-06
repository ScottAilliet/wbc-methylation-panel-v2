# Materials and Methods

## 1. Data Source

Whole-genome bisulfite sequencing (WGBS) data from the Loyfer et al. human methylation atlas
(GEO accession GSE186458) [1]. This atlas comprises 207 samples spanning 82 human cell types,
including 7 immune cell types of interest (monocytes, B cells, NK cells, granulocytes, CD3 T cells,
CD8 T cells, CD4 T cells) and 75 non-immune cell types (liver, heart, kidney, brain, pancreas,
lung, colon, breast, skin, muscle, etc.).

Reference genome: hg19 (GRCh37).

## 2. DMR Discovery

### 2.1. Tools

Differentially methylated regions (DMRs) were identified using `wgbstools find_markers`
(https://github.com/nloyfer/wgbs_tools), which compares methylation levels between target
and background cell types across pre-defined genomic blocks.

### 2.2. Block Definition

Blocks were defined using the GSE186458_blocks.s205.bed.gz file (205 CpG blocks), provided
with the atlas. Each block represents a contiguous genomic region with a defined set of CpG sites.

### 2.3. Target and Background Definition

For each immune cell type, DMRs were identified as blocks that are:
- **Hypomethylated** in the target cell type (target mean methylation < 0.25 quantile)
- **Methylated** in all background cell types (background mean methylation > 0.975 quantile)
- Differentially methylated with delta_means ≥ 0.2 (non-T cell types) or ≥ 0.1 (T cell subtypes)

Background definitions:
- **Monocytes, B cells, NK, Granulocytes**: All 81 other human cell types
- **CD3 T cells (pan-T)**: All 73 non-T cell types
- **CD8 T cells**: Non-T + CD3 + CD4 lineage (78 groups)
- **CD4 T cells**: Non-T + CD3 + CD8 lineage (78 groups)

### 2.4. Filtering Parameters

| Parameter | Value |
|---|---|
| Minimum CpGs per block | 5 |
| Preferred minimum CpGs | 7 |
| Minimum coverage | 5 |
| P-value threshold | 0.05 |
| NA rate (target) | 0.334 |
| NA rate (background) | 0.334 |
| Target quantile | 0.25 |
| Background quantile | 0.025 |
| Top markers per cell type | 300 |
| Sort by | delta_means |
| Direction | Hypomethylated only |

### 2.5. Per-CpG Methylation Extraction

For each DMR block, individual CpG positions and methylation values were extracted:
- CpG genomic positions from the hg19 CpG index (CpG.bed.gz)
- Per-CpG beta values from each of the 207 beta files (binary format: 2 bytes per CpG,
  uint8 methylated count + uint8 total count)
- Each CpG labeled (CpG_1, CpG_2, ...) for primer design linkage
- Target mean beta and background mean beta calculated per CpG

## 3. Primer Design

### 3.1. Bisulfite Conversion

Genomic sequences were bisulfite-converted in silico to produce 6 strands:
- **S**: Sense (original top strand)
- **SM**: Sense methylated (C→T except CpG C's preserved)
- **SU**: Sense unmethylated (all C→T)
- **A**: Antisense (original bottom strand, reverse complement)
- **AM**: Antisense methylated
- **AU**: Antisense unmethylated

Methylation-specific PCR (MSP) primers were designed on the appropriate converted strand:
- Methylated primers: from SM or AM
- Unmethylated primers: from SU or AU

### 3.2. Primer3 Configuration

Primers were designed using primer3-py v2.3.0 with Primer3Plus-compatible settings
(161 parameters). Key parameters:

| Parameter | Value |
|---|---|
| Tm range | 58–62°C (user-specified, required) |
| Primer size | 18–23 bp (optimal: 20) |
| Product size | 60–150 bp |
| GC content | 20–80% |
| Max self-complementarity | 8.00 |
| Max self-end complementarity | 3.00 |
| Max hairpin Tm | 47°C |
| Max pair Tm difference | 1°C |
| GC clamp | 1 |
| Max end GC | 2 |
| Max end stability | 9.0 kcal/mol |
| Max poly-X | 5 |
| Salt correction | SantaLucia (method 1) |
| Tm formula | SantaLucia 1998 (method 1) |

### 3.3. Bisulfite-Specific Constraints

- Minimum 2 CpGs per primer
- Minimum 4 CpGs per primer pair
- Terminal CpG at 3' end preferred for methylation discrimination

### 3.4. Salt Conditions

Two salt presets are available:
- **Primer3Plus**: 50mM Na+, 3mM Mg2+, 1.2mM dNTP, 250nM primer
- **Roche DLC**: 100mM Na+, 4.5mM Mg2+, 1.2mM dNTP, 250nM primer

Tm values are recalculated according to the selected salt conditions.

## 4. QC Filters

### 4.1. Bowtie2 Specificity (6-Strand)

Each primer pair was screened against 6 genome states:
1. Unconverted genome (sense + antisense)
2. Converted unmethylated genome (sense + antisense)
3. Converted methylated genome (sense + antisense)

A primer pair passes if it aligns uniquely to the intended genome state with no
significant off-target alignments (≤4 mismatches) in other states.

### 4.2. Secondary Structure Screening

Hairpin formation MFE was calculated using primer3.calc_hairpin with SantaLucia
thermodynamic parameters [2]:
- mv_conc = 50 mM, dv_conc = 1.5 mM, dntp_conc = 0.6 mM, dna_conc = 50 nM, temp = 37°C
- Cutoff: MFE ≥ -1.5 kcal/mol (conservative)

### 4.3. Common SNP Screening

Primers were screened against dbSNP common variants (MAF ≥ 1%):
- 2-digit score: first digit = SNP in primer body, second digit = SNP in 3' end (last 5 nt)
- Score 00 = clean, 10 = body SNP, 01 = 3' SNP, 11 = both

### 4.4. Primer-Dimer Prediction (DimerDetective)

Primer-dimer formation risk was assessed using DimerDetective, a thermodynamic scoring
method based on nearest-neighbour 3'-end stability calculated with primer3.calc_end_stability
(SantaLucia 1998) [3]. For each primer pair, both heterodimer orientations are evaluated:
the 3' end of the left primer annealed to the right primer, and vice versa. The minimum ΔG
across both orientations (end_min_dg, kcal/mol) determines the risk classification.

DimerDetective was validated on 200 empirically validated qPCR assays (100 dimer-forming,
100 clean) and applies three risk tiers:
- **High**: end_min_dg ≤ −2.48 kcal/mol (100% specificity)
- **Medium**: −2.48 < end_min_dg ≤ −0.18 kcal/mol (mixed zone)
- **Low**: end_min_dg > −0.18 kcal/mol (100% sensitivity)

Conservative cutoff: −1.0 kcal/mol (95% sensitivity, 55% specificity).

Thermodynamic parameters: mv_conc = 50 mM, dv_conc = 1.5 mM, dntp_conc = 0.6 mM,
dna_conc = 50 nM, temp_c = 37°C.

**Caveat**: DimerDetective was validated on standard qPCR primers with 40–60% GC content.
Bisulfite-converted MSP primers exhibit substantially lower GC content (typically 10–25%),
which falls outside the validation range. Risk labels and end_min_dg values for MSP primers
should be interpreted as indicative rather than quantitative.

## 5. CFF Assays

C-free fragment (CFF) assays were designed for gDNA extraction QC (fragmentation and
recovery measurement). C-free regions (≥120 bp, no C on the top strand) were identified
by scanning the hg19 genome. Primers were designed for two amplicon sizes:
- **CFF_50**: 50 bp amplicon
- **CFF_100**: 100 bp amplicon

FWD primers contain only A/T/G (no C), REV primers contain only A/T/C (no G).

## 6. Conversion Controls

Housekeeping gene regions were screened for always-unmethylated status across all 207
atlas samples. Candidate genes: GAPDH, RPP30, HPRT1, B2M, ACTB. Selection criteria:
- Mean methylation < 0.1 across all cell types
- No pseudogenes
- Consistent unmethylation in all samples

## 7. Output Format

### 7.1. XLSX Output

Excel file with 27 columns per primer pair, matching the U-assays reference format:
assay_id, seq_id, template_used, assay, left/right_primer_display, Tm, product_size,
CpG counts, mismatch scores, bowtie_passes_filter, structure MFE, dimer prediction,
common_variant_score, mapping_error_note.

### 7.2. PDF Output

One page per primer pair with:
- 6-strand sequence visualization (S, SM, SU, A, AM, AU)
- Primer information table
- Assay characteristics
- QC filter results
- DimerDetective calculations footer

## References

1. Loyfer, N. et al. A human methylation atlas for normal and cancer cell types. Nature 617, 690–697 (2023).
2. SantaLucia, J. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proc. Natl. Acad. Sci. USA 95, 1460–1465 (1998).
3. SantaLucia, J. & Hicks, D. The thermodynamics of DNA structural motifs. Annu. Rev. Biophys. Biomol. Struct. 33, 415–440 (2004).
