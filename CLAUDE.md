# CLAUDE.md — BAMVIS-GENE

## Project Overview

BAMVIS-GENE is a bioinformatics visualization tool that renders RNA-seq BAM alignment data as SVG images at the gene level. It draws gene models (CDS, UTR, mRNA) from GFF3 annotations and overlays aligned sequencing reads from BAM files, showing splice junctions, paired-end bridges, and strand orientation.

The tool targets plant genomics data (Arabidopsis thaliana, Chlamydomonas reinhardtii, Eutrema salsugineum) using Phytozome GFF3 annotation files.

## Quick Start

```sh
# Install dependencies
pip install -r requirements.txt

# Step 1: Index a GFF3 file
python bamvis.py index Athaliana_167_TAIR10.gene.gff3

# Step 2: Draw a gene
python bamvis.py draw --pickle Athaliana_167_TAIR10.gene.gff3.pandas.df.pk \
                      --bam sample.sorted.bam \
                      --gene AT1G01010

# Step 2: Draw a genomic region (no gene model)
python bamvis.py draw --pickle Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk \
                      --bam reads.sorted.bam \
                      --region Chr9:3586521-3608027 --strand +
```

## Repository Structure

```
BAMVIS-GENE/
├── bamvis.py                                            # Unified CLI tool (Python 3, argparse)
├── requirements.txt                                     # pip dependencies
├── 1.RNAseq.bam.visualization.gffindexing.py          # Legacy Step 1: GFF3 → pickle (longest transcript filter)
├── 1.RNAseq.bam.visualization.gffindexing.nonlongest.py # Legacy Step 1 variant: no longest transcript filter
├── 2.RNAseq.bam.visualization.draw.py                  # Legacy Step 2: samtools-based SVG drawer
├── 2.RNAseq.bam.visualization.draw_old.py              # Legacy Step 2: older copy
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.py           # Legacy: pysam, single gene, hardcoded paths
├── 2.RNAseq.bam.visualization.draw.pysam.transcript.underconstruction.py # Legacy: pysam, transcript-level indexing
├── 2.RNAseq.bam.visualization.draw.pysam.transcript.underconstruction.all.py # Legacy: pysam batch mode
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.region.py    # Legacy: pysam, arbitrary region
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.region.notebook.ipynb # Jupyter notebook (region mode)
├── kang.py                                              # Shared utility module (bioinformatics helpers)
├── *.gff3                                               # GFF3 annotation files (large, Phytozome format)
├── *.pandas.df.pk                                       # Pre-built pandas pickle index files
├── *.svg / *.png                                        # Example output visualizations
├── temp.sam / temp.sam.cut                              # Temporary intermediate SAM files
└── README.md
```

## Unified CLI: bamvis.py

`bamvis.py` consolidates all script variants into a single Python 3-compatible CLI with two subcommands.

### `bamvis.py index` — GFF3 Indexing

```sh
python bamvis.py index <gff3_file> [--no-longest]
# Output: <gff3_file>.pandas.df.pk
```

Options:
- `--no-longest` — Keep all transcripts instead of filtering for the longest per gene

### `bamvis.py draw` — BAM Visualization

```sh
# Gene mode (draws gene model + reads)
python bamvis.py draw --pickle <pk> --bam <bam1> [--bam <bam2>] --gene <gene_name>

# Region mode (draws reads only, no gene model)
python bamvis.py draw --pickle <pk> --bam <bam1> --region chr:start-end [--strand +/-]
```

Options:
- `--pickle` — Pandas pickle index file (`.pandas.df.pk`)
- `--bam` — BAM file path (repeat for multiple BAM files)
- `--gene` — Gene name to visualize
- `--region` — Genomic region as `chr:start-end`
- `--strand` — Strand for region mode (`+` or `-`, default: `+`)
- `--rows` — Max read packing rows (default: 20)
- `--output` / `-o` — Output SVG filename (auto-generated if omitted)

## Legacy Scripts

The original numbered scripts (`1.*.py`, `2.*.py`) remain in the repository. They are Python 2 scripts with hardcoded paths. See the sections below for details if you need to work with them.

### Legacy Step 1 — GFF Indexing (`1.*.py`)

```sh
python 1.RNAseq.bam.visualization.gffindexing.py <gff3_file>
```

Two variants:
- **gffindexing.py** — Indexes by `(genename, longest)`, filtering for longest transcript per gene
- **gffindexing.nonlongest.py** — Indexes by `(genename)` only, keeping all transcripts

### Legacy Step 2 — BAM Visualization (`2.*.py`)

**Original (samtools-based):**
```sh
python 2.RNAseq.bam.visualization.draw.py <pickle_file> <bam_file> <gene_name>
```

**pysam-based variants** have hardcoded reference paths and BAM lists inside the scripts. Edit the script globals before running.

## Dependencies

```sh
pip install -r requirements.txt
```

- `pandas` — DataFrame manipulation and GFF3 parsing
- `pysam` — BAM file reading
- `numpy` — 2D array for read layout/packing
- `tqdm` — Progress bars during BAM iteration

Legacy scripts additionally require:
- **Python 2** with `from __future__ import print_function`
- `samtools` on PATH (for the original non-pysam drawer)
- `kang.py` importable (some scripts use `sys.path.append`)

## Utility Module: kang.py

Provides bioinformatics helper functions:
- `get_block(array, depth_cut, lim_len_block)` — Find contiguous blocks above a depth threshold
- `flagparser(intFlag)` — Parse SAM bitwise flags into a human-readable dict
- `rev_comp(strSeq)` — Reverse complement a DNA sequence
- `translation(strSeq)` — Translate DNA to protein using the standard genetic code
- `Fasta2dic(file_fasta)` / `fasta2dic(file_fasta, dic)` — Parse FASTA files into dicts
- `dic2fa(dic, filename)` — Write a dict back to FASTA format
- `list2txt(outfilename, inlist)` — Write a list to a text file

## SVG Output Conventions

- **Canvas**: 1000px wide; height scales with number of reads and BAM files
- **Gene model**: Drawn at the top — mRNA as a thin line, CDS as colored rectangles, UTR as colored rectangles
- **Reads**: Packed vertically using a greedy row-packing algorithm (`gene_space` 2D array)
- **Paired-end bridge**: V-shaped line connecting mate pairs
- **Splice junctions**: Thin lines for N (intron skip) and I (insertion) CIGAR operations
- **Color scheme** ("Haruhi theme"):
  - UTR: `#f50002` (red)
  - CDS: `#53a1b5` (teal)
  - Forward strand reads: `#5d3c2d` (brown)
  - Reverse strand reads: `#f5b024` (gold)
  - Legacy samtools drawer uses: CDS=yellow, UTR=blue, reads=pink

## Key Conventions for AI Assistants

- **bamvis.py is the primary entry point**: Use `bamvis.py` for new work. It is Python 3 compatible and uses `argparse` with no hardcoded paths.
- **Legacy scripts are Python 2**: The `1.*.py` and `2.*.py` files use `from __future__ import print_function`, old pandas API (`df.sort()` instead of `df.sort_values()`, `df.sortlevel()` instead of `df.sort_index()`), and `dict.keys()` returning a list.
- **No test suite or CI/CD**: There are no automated tests. Verify changes manually.
- **Hardcoded paths in legacy scripts**: Many pysam-based scripts contain hardcoded absolute paths (e.g., `/ref/analysis/...`, `/mnt/c/ubuntu.download/...`). These are not present in `bamvis.py`.
- **Script naming**: Legacy files are numbered by pipeline step (`1.` = indexing, `2.` = drawing). Variants are appended as dot-separated descriptors.
- **Large binary files in repo**: GFF3 files, pickle files, SAM files, and images are committed directly. Be mindful of repo size.
- **Coordinate system**: GFF3 uses 1-based coordinates. A 100bp padding is added on each side of the gene region for visualization.
- **CIGAR parsing**: Uses regex `r'(\d+)(\w)'` to parse CIGAR strings. Handles M (match), N (intron skip), and I (insertion) operations.
- **Read filtering**: Skips improperly paired reads, duplicates, QC-failed reads, and secondary alignments.
