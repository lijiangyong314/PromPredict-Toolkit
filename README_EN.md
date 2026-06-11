# PromPredict Toolkit

> **Physics over black boxes.**  
> Promoter prediction based on DNA duplex thermodynamic stability.

[![Python](https://img.shields.io/badge/python-3.7%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2Fb906535k-red)](https://doi.org/10.1039/B906535K)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20wsl-lightgrey)]()

**GitHub:** [github.com/lijiangyong314/PromPredict-Toolkit](https://github.com/lijiangyong314/PromPredict-Toolkit)

---

## Overview

**PromPredict** predicts promoter regions in genomic DNA based on a simple yet powerful principle: **promoter DNA must be locally unstable to allow RNA polymerase to melt the double helix and initiate transcription.**

Unlike machine-learning "black boxes," PromPredict computes the average free energy (ΔG) of every 15-bp sliding window along the input sequence using the Breslauer dinucleotide parameters. A region where the ΔG rises significantly above the flanking baseline is flagged as a candidate promoter — because nature makes promoters easier to melt.

- **No training data** — thresholds are derived from experimentally validated TSS datasets in *E. coli*, *B. subtilis*, and *M. tuberculosis*
- **No GPU** — pure physics, runs on a laptop CPU
- **No model files** — 2.4 MB standalone binary
- **Interpretable** — every prediction comes with a Dmax (stability difference) you can reason about

This toolkit wraps the official IISc PromPredict binary with:
- A production-grade Python visualization suite (`plot_prompredict.py`)
- A comprehensive command-line usage guide
- QC scripts and confidence-level filtering

| Method | Approach | ML Model | GPU Required | Interpretable |
|--------|----------|----------|-------------|---------------|
| **PromPredict** | Thermodynamics (ΔG) | None | No | ✅ Yes (Dmax) |
| iPro-MP | BERT Transformer | 36.5 GB | Yes | ❌ No |
| Promotech | Random Forest / RNN | Custom | Optional | Partially |

---

## Features

- **Zero-dependency binary**: download and run — no Python, no R, no Conda required
- **Multisequence mode**: batch process hundreds of genes in one run
- **Genome mode**: scan entire bacterial chromosomes (>10 MB) in sliding-window mode
- **GC-adaptive thresholds**: automatically adjusts free energy cutoffs based on GC content
- **Confidence scoring**: every prediction comes with a `Dmax` value (higher = stronger signal)
- **Publication-ready visualization**: `plot_prompredict.py` generates one profile plot per sequence with promoter regions highlighted (multisequence mode). Genome-mode browser-track visualization is planned
- **QC pipeline**: built-in scripts to check coverage, Dmax distribution, and format integrity (both modes)

---

## Quick Start

```bash
# 1. Download
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_mulseq" \
  -O prompredict && chmod +x prompredict

# 2. Prepare input (FASTA)
echo -e ">gene1\nATCGATCGATCG..." > test.fasta

# 3. Run
printf "test.fasta\n100\ndefault\n" | ./prompredict

# 4. Visualize
python3 plot_prompredict.py test
```

That's it. Three output files are generated: a promoter list, a free-energy profile, and GC statistics. The visualization script produces one high-resolution PNG per sequence plus a confidence summary chart.

**Genome-mode variant:**

```bash
# Download the genome binary
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_genome_V2" \
  -O prompredict_genome && chmod +x prompredict_genome

# Run on a bacterial chromosome
printf "genome.fna\n100\ndefault\n" | ./prompredict_genome
```

> ⚠️ Genome-mode output has a **different format** (11 columns per record vs. 8 in multisequence mode). `plot_prompredict.py` currently supports multisequence mode only. See [Genome-Mode QC](#genome-mode-qc) for manual filtering commands.

---

## Installation

### Prerequisites

- **OS**: Linux (x86-64) or WSL2 on Windows
- **Disk**: < 5 MB for the binary; ~1 MB per 100 sequences for output
- **Python** (optional, for visualization): 3.7+ with `matplotlib` and `numpy`

### Step 1 — Download the binary

```bash
mkdir -p ~/bio_tools/prompredict && cd ~/bio_tools/prompredict

# Multisequence mode (most common)
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_mulseq" \
  -O prompredict
chmod +x prompredict

# Optional: whole-genome mode (for >10 MB FASTA)
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_genome_V2" \
  -O prompredict_genome
chmod +x prompredict_genome
```

> **Note**: The IISc server uses an expired SSL certificate (RapidSSL, issued 2011). The `--no-check-certificate` flag is **required**.

### Step 2 — Install Python visualization dependencies (optional)

```bash
# Using your existing environment
pip install matplotlib numpy

# Or create an isolated environment
conda create -n prompredict python=3.9 matplotlib numpy -y
conda activate prompredict
```

### Step 3 — Verify

```bash
./prompredict 2>&1 | head -3
# Expected output:
# PromPredict : Date of Creation : 20 March 2008
# Predicts promoter region over given set of sequences based on relative stability
```

---

## Input Format

PromPredict accepts **multi-FASTA** files. Each entry must have a header line starting with `>` followed by a sequence identifier, and the DNA sequence on subsequent lines.

```fasta
>gene1_upstream_500bp
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG...
>gene2_upstream_500bp
GGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC...
```

| Requirement | Details |
|-------------|---------|
| **Minimum length** | 200 bp (recommended: 300–500 bp) |
| **Maximum sequences** | No hard limit (tested to 10,000+) |
| **Allowed characters** | A, T, G, C (uppercase recommended; case-insensitive in practice) |
| **Gaps/ambiguous bases** | Not supported — remove or replace before running |

---

## Output Format

After running, three files are generated. Replace `<prefix>` with your input filename (without `.fasta`).

### 1. `<prefix>_PPde.txt` — Predicted promoter regions

The primary output. Each predicted promoter spans **two lines**:

```
ID    <seq_id>    <GC_content>
><start>..<end>    <length>    <sequence>    <lsp>    <lspe>    <dmax_pos>    <Dmax>    <Dave>
```

**Example** (the *E. coli* ompF promoter):

```
ID    PM00066ompFpTU00050986315-    38.76
>301..402    102    ttaatttaaatataaggaaatcat...    308    -11.96    328    2.34    1.86
```

| Column | Field | Description |
|--------|-------|-------------|
| 1 | `>start..end` | Promoter coordinates (**1-based**) |
| 2 | `length` | Promoter region length (bp) |
| 3 | `sequence` | DNA sequence of predicted promoter |
| 4 | `lsp` | **Local Stability Peak** — position of maximum instability; approximates the −10 element |
| 5 | `lspe` | Energy at LSP (kcal/mol) |
| 6 | `dmax_pos` | Position of maximum stability difference |
| 7 | **`Dmax`** | **Stability difference** (promoter vs. flanking). **Primary confidence metric.** |
| 8 | `Dave` | Average stability difference across the promoter |

#### Dmax Confidence Guide

| Dmax | Confidence | Recommended Action |
|------|-----------|-------------------|
| **≥ 3.5** | 🔥 Very High | Experimental validation |
| **3.0 – 3.5** | ✅ High | Prioritize for validation |
| **2.0 – 3.0** | ⚠️ Medium | Cross-validate with iPro-MP or motif search |
| **1.5 – 2.0** | ❓ Weak | Requires additional evidence |
| **< 1.5** | ❌ Very Low | Likely noise; ignore |

> 💡 **Coordinate system**: PromPredict uses **1-based indexing** (position 1 = first nucleotide). When extracting substrings in Python (`seq[start:end]`), convert with `start_py = start - 1`.

### 2. `<prefix>_stb.txt` — Free-energy profile

Each line represents one 15-bp sliding window:

```
ID    <seq_id>    ΔG₁    ΔG₂    ΔG₃    ...
```

- **More negative** (e.g., −22 kcal/mol): very stable duplex → unlikely promoter
- **Less negative** (e.g., −12 kcal/mol): unstable duplex → candidate promoter region

The number of rows is `sequence_length − 14` (because 15-bp windows lose 7 bp of coverage at each end).

### 3. `<prefix>_GCstat.txt` — GC-content distribution

```
35-40% GC: 9 sequences
45-50% GC: 1 sequence
```

Used internally for threshold selection. You generally do not need to inspect this file.

### Genome-Mode Output (Format Differences)

When using `PromPredict_genome` instead of `PromPredict_mulseq`, the output file structure is modified:

**`_PPde.txt` — flat 11-column format (no `ID`/`>` headers)**

```
# Start of the 1000nt window  %GC  pstart  pend  length  promoter_seq  lsp  lspe  Dmax_pos  Dmax  Dave
0    15.90   500   607   108   tcgcaaat...   549   -11.42   555   1.85   1.52
750  15.40   1173  1295  123   tttttataa...  1244  -11.84   1239  1.46   1.16
```

| Column | Index | Multisequence equivalent |
|--------|-------|--------------------------|
| Dmax | **9** (not 6) | `parts[6]` |
| Dave | **10** (not 7) | `parts[7]` |
| `pstart` / `pend` | 2 / 3 | `>start..end` (col 0) |

**`_stb.txt` — uses window start position instead of sequence ID**

```
# Start of the 1000nt window  ΔG₁  ΔG₂  ΔG₃  ...
0    -13.73  -13.85  -13.70  ...
750  -14.25  -14.53  -14.13  ...
```

> 🔧 If you need genome-mode visualization, use the QC commands in [Genome-Mode QC](#genome-mode-qc) to extract coordinates manually, or open an issue to track genome-mode plot support.

---

## Usage

### Interactive mode

PromPredict is interactive. It prompts for three parameters in sequence:

```
Enter the Input File Name:        your_file.fasta
Enter the E1 region window size:  100
Whole genome GC content:          default
```

### Non-interactive mode (recommended for scripts)

```bash
printf "input.fasta\n100\ndefault\n" | ./prompredict
```

The `printf` command feeds three lines into the interactive prompts via a Unix pipe.

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Input file** | *(required)* | Multi-FASTA file |
| **E1 window** | `100` | Upstream comparison window (nt). `100` = high sensitivity + precision. Use `50` for second-pass if `100` misses known promoters |
| **GC content** | `default` | `default` = auto-detect from 1000-nt fragments. Alternatively, specify the known whole-genome GC% (e.g., `50.5`) |

### Whole-genome scanning

For genomes > 10 MB, use the genome-mode binary:

```bash
printf "genome.fna\n100\ndefault\n" | ./prompredict_genome
```

Two versions are available:
- `PromPredict_genome_V1` — for genomes < 10 MB
- `PromPredict_genome_V2` — for genomes ≥ 10 MB

---

## Quality Control

Always run these QC checks after prediction:

```bash
# 1. Coverage: did all sequences get analyzed?
grep -c "^>" input.fasta                    # Input sequences
grep -c "^ID" input_PPde.txt                # Sequences with ≥1 prediction

# 2. Dmax distribution: are values in a healthy range?
awk -F'\t' '/^>/ {print $7}' input_PPde.txt | sort -nr | head -20

# 3. Anomalies: impossible Dmax values?
awk -F'\t' '/^>/ && ($7 > 10 || $7 < 0) {print}' input_PPde.txt

# 4. Breakdown by confidence tier
high=$(awk -F'\t' '/^>/ && $7>=3.0' input_PPde.txt | wc -l)
med=$(awk  -F'\t' '/^>/ && $7>=2.0 && $7<3.0' input_PPde.txt | wc -l)
low=$(awk  -F'\t' '/^>/ && $7<2.0' input_PPde.txt | wc -l)
echo "High: $high | Medium: $med | Low: $low"
```

If **all** sequences show Dmax < 2.0:
1. Check input sequence length (should be ≥ 200 bp)
2. For high-GC genomes (>60%), manually specify the GC% value
3. Try re-running with `E1=50` for shorter promoter detection

### Genome-Mode QC

Genome-mode `_PPde.txt` has a **flat format** (no `ID`/`>` lines) with **11 tab-separated columns**:

```
col 0: window_start  col 3: pstart  col 4: pend  col 9: Dmax  col 10: Dave
```

```bash
# Genome-mode QC commands (different column indices!)
awk -F'\t' 'NR>2 && $10 >= 3.0' genome_PPde.txt | wc -l     # High confidence
awk -F'\t' 'NR>2 && $10 >= 2.0 && $10 < 3.0' genome_PPde.txt | wc -l  # Medium
awk -F'\t' 'NR>2 && $10 < 2.0' genome_PPde.txt | wc -l      # Low

# Extract all high-confidence predictions with coordinates
awk -F'\t' 'NR>2 && $10 >= 3.0 {print $1"\t"$4"\t"$5"\t"$9"\t"$10}' genome_PPde.txt > high_conf.bed
```

---

## Visualization

The included `plot_prompredict.py` script generates publication-quality figures:

```bash
python3 plot_prompredict.py <prefix>
```

> ⚠️ **Currently supports multisequence mode only.** Genome-mode visualization (genome-browser-style ΔG tracks) is planned for a future release. For genome-mode results, use the [Genome-Mode QC](#genome-mode-qc) commands to filter and export high-confidence predictions.

**Output directory** (`<prefix>_profiles/`):

| File | Description |
|------|-------------|
| `<prefix>_confidence_summary.png` | Horizontal bar chart: all predictions ranked by Dmax, color-coded by confidence |
| `<seq_id_1>.png` | Free-energy profile for sequence 1, with predicted promoters highlighted |
| `<seq_id_2>.png` | ... one PNG per sequence |

Each profile plot shows:
- Blue line: 15-bp sliding-window ΔG curve
- Red dashed line: sequence-wide mean ΔG
- Colored bands: predicted promoter regions (green = high confidence, orange = medium, red = low)
- Dmax annotations above each band

### Example

```bash
# After running PromPredict on sample data:
printf "sample_seq.fasta\n100\ndefault\n" | ./prompredict
python3 plot_prompredict.py sample_seq

# Generates:
#   sample_seq_profiles/
#   ├── sample_seq_confidence_summary.png   (overview)
#   ├── PM0-4641yaiS...png                  (per-sequence profiles)
#   └── ...
```

---

## Parameter Tuning

| Symptom | Biological Explanation | Fix |
|---------|----------------------|-----|
| **Too many false positives** (5+ promoters per gene) | Threshold is too low; ordinary unstable regions flagged as promoters | Raise Dmax cutoff to **3.5** during filtering |
| **Too many false negatives** (known promoters missed) | E1=100 is too broad, masking short promoter signals | Re-run with **E1=50** |
| **High-GC genome** (GC > 60%) | Global ΔG baseline is much lower; signal-to-noise ratio degrades | Specify **actual GC%** (not `default`) and relax Dmax cutoff to **>2.0** |
| **Low-GC genome** (GC < 30%) | Global ΔG baseline is higher; may miss weaker signals | Specify **actual GC%** |

---

## Integration with AI Tools

PromPredict excels at **explainable first-pass screening**. For a complete promoter annotation pipeline, combine it with complementary tools:

| Stage | Tool | Why It's Irreplaceable |
|-------|------|----------------------|
| **1. Screen** | **PromPredict** | Only physics-based method. Every prediction is explainable: "DNA here is easier to melt." |
| **2. Validate** | iPro-MP / DeePromoter | BERT/CNN models recognize complex σ-factor-specific sequence motifs that PromPredict cannot distinguish |
| **3. Score** | promoterCalculator | Quantifies −35/−10 element match against position weight matrices |
| **4. Terminate** | cryptkeeper / transtermHP | Rho-dependent and independent terminator prediction — often co-localized with promoters |
| **5. Structure** | ViennaRNA | RNA secondary structure at the promoter can affect transcription initiation efficiency |

---

## Scope and Limitations

### Best Performance
- Gram-negative bacteria (especially *Proteobacteria*) with GC 35–55%
- Upstream regions of protein-coding genes (200–500 bp)
- σ70-type promoters (the predominant type in *E. coli*)

### Use with Caution
- High-GC bacteria (*Streptomyces* >72%, *Mycobacterium* >65%) or Archaea
- σ54-type promoters (different melting mechanism; lower sensitivity)
- GC < 30% genomes (may miss weak signals)

### Not Recommended
- Mammalian genomes — CpG islands, chromatin structure, and distal enhancers are beyond PromPredict's design
- Predicting promoter **strength** — PromPredict predicts **location**, not expression level. Use `promoterCalculator` for strength

---

## FAQ

<details>
<summary><b>Why is the _stb.txt line count less than the sequence length?</b></summary>

_stb.txt records values for 15-bp sliding windows. Row N corresponds to the **center** of window N, so there are `sequence_length − 14` rows total (each end loses 7 bp of coverage).
</details>

<details>
<summary><b>The predicted promoter overlaps a CDS. Is this an error?</b></summary>

Not necessarily. PromPredict identifies **thermodynamically unstable regions**, which can occur inside genes: internal promoters, Rho-dependent terminator regions, or horizontal gene transfer insertion sites. Cross-validate with −10/−35 box motif search (<code>TATAAT</code> / <code>TTGACA</code>).
</details>

<details>
<summary><b>How do I choose between multisequence and genome mode?</b></summary>

- **Multisequence** (`PromPredict_mulseq`): up to a few hundred gene upstream regions, each in a separate FASTA entry. Ideal for targeted analysis.
- **Genome** (`PromPredict_genome_V2`): one large FASTA (complete chromosome or scaffold). Ideal for whole-genome promoter scanning.
</details>

<details>
<summary><b>PromPredict reports a promoter on every sequence. Is that normal?</b></summary>

It can happen with sequences longer than 500 bp — the binary applies a permissive first pass. Filter by Dmax: keep only Dmax ≥ 2.0 (or ≥ 3.0 for high-stringency analysis).
</details>

<details>
<summary><b>Can I run this on macOS?</b></summary>

The IISc binary is compiled for Linux x86-64. On macOS, you have two options:
1. **Docker**: `docker run -v $(pwd):/data ubuntu:20.04 /data/prompredict ...`
2. **Web server**: use the online version at https://dna.mbu.iisc.ac.in/prompredict/
</details>

---

## Citation

If you use PromPredict in your research, please cite:

> Rangannan, V., & Bansal, M. (2009).  
> **Relative stability of DNA as a generic criterion for promoter prediction: whole genome annotation of microbial genomes with varying nucleotide base composition.**  
> *Molecular BioSystems*, 5(12), 1758–1769.  
> DOI: [10.1039/B906535K](https://doi.org/10.1039/B906535K)

For eukaryotic applications, also cite:

> Rangannan, V., & Bansal, M. (2018).  
> **Identification of putative promoters in 48 eukaryotic genomes on the basis of DNA free energy.**  
> *Scientific Reports*, 8, 4520.  
> DOI: [10.1038/s41598-018-22129-8](https://doi.org/10.1038/s41598-018-22129-8)

---

## Contributing

Bug reports, feature requests, and pull requests are welcome at [github.com/lijiangyong314/PromPredict-Toolkit](https://github.com/lijiangyong314/PromPredict-Toolkit). Please open an issue before submitting a PR for major changes.

## License

This toolkit (Python scripts, documentation, visualization code) is released under the **MIT License**. See [LICENSE](LICENSE) for details.

The PromPredict binary is © Molecular Biophysics Unit, Indian Institute of Science. Redistribution terms are available at the [official download page](https://dna.mbu.iisc.ac.in/prompredict/download.html).

## Acknowledgments

- **Vetriselvi Rangannan & Manju Bansal** (IISc) — for developing PromPredict and making it freely available
- The Breslauer dinucleotide free-energy parameters (Breslauer et al., 1986 *PNAS*)
- The experimentally validated TSS datasets from *E. coli* K-12, *B. subtilis*, and *M. tuberculosis* used for threshold derivation
