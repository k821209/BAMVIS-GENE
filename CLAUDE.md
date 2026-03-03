# CLAUDE.md — BAMVIS-GENE

## Project Overview

BAMVIS-GENE is a bioinformatics visualization tool that renders RNA-seq BAM alignment data as SVG images at the gene level. It draws gene models (CDS, UTR, mRNA) from GFF3 annotations and overlays aligned sequencing reads from BAM files, showing splice junctions, paired-end bridges, and strand orientation.

The tool targets plant genomics data (Arabidopsis thaliana, Chlamydomonas reinhardtii, Eutrema salsugineum) using Phytozome GFF3 annotation files.

## Repository Structure

```
BAMVIS-GENE/
├── 1.RNAseq.bam.visualization.gffindexing.py          # Step 1: GFF3 → pandas pickle (longest transcript filter)
├── 1.RNAseq.bam.visualization.gffindexing.nonlongest.py # Step 1 variant: no longest transcript filter
├── 2.RNAseq.bam.visualization.draw.py                  # Step 2: Original samtools-based SVG drawer
├── 2.RNAseq.bam.visualization.draw_old.py              # Step 2: Older copy of samtools-based drawer
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.py           # Step 2: pysam-based, single gene, hardcoded paths
├── 2.RNAseq.bam.visualization.draw.pysam.transcript.underconstruction.py # Step 2: pysam-based, transcript-level indexing
├── 2.RNAseq.bam.visualization.draw.pysam.transcript.underconstruction.all.py # Step 2: pysam batch mode (all transcripts)
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.region.py    # Step 2: pysam-based, arbitrary genomic region
├── 2.RNAseq.bam.visualization.draw.pysam.underconstruction.region.notebook.ipynb # Jupyter notebook version (region mode)
├── kang.py                                              # Shared utility module (bioinformatics helpers)
├── *.gff3                                               # GFF3 annotation files (large, Phytozome format)
├── *.pandas.df.pk                                       # Pre-built pandas pickle index files
├── *.svg / *.png                                        # Example output visualizations
├── temp.sam / temp.sam.cut                              # Temporary intermediate SAM files
└── README.md
```

## Two-Step Pipeline

### Step 1 — GFF Indexing (`1.*.py`)

Parses a GFF3 file into a pandas DataFrame indexed by gene name, serialized as a pickle (`.pandas.df.pk`).

```sh
python 1.RNAseq.bam.visualization.gffindexing.py <gff3_file>
# Output: <gff3_file>.pandas.df.pk
```

Two variants:
- **gffindexing.py** — Indexes by `(genename, longest)`, filtering for the longest transcript per gene
- **gffindexing.nonlongest.py** — Indexes by `(genename)` only, keeping all transcripts

### Step 2 — BAM Visualization (`2.*.py`)

Reads the pickle index + BAM file(s), renders SVG output with gene model and aligned reads.

**Original (samtools-based):**
```sh
python 2.RNAseq.bam.visualization.draw.py <pickle_file> <bam_file> <gene_name>
```

**pysam-based variants** have hardcoded reference paths and BAM lists inside the scripts. Edit the script globals before running:
```sh
python 2.RNAseq.bam.visualization.draw.pysam.underconstruction.py <gene_name>
```

**Region-based variant:**
```sh
python 2.RNAseq.bam.visualization.draw.pysam.underconstruction.region.py <chromosome> <start> <end> <strand>
```

## Dependencies

- **Python 2** (uses `from __future__ import print_function`, older pandas API like `df.sort()`)
- `pandas` — DataFrame manipulation and GFF3 parsing
- `pysam` — BAM file reading (install: `pip2 install pysam`)
- `numpy` — 2D array for read layout/packing
- `tqdm` — Progress bars during BAM iteration
- `samtools` — Required on PATH for the original (non-pysam) drawer
- `kang.py` — Local utility module (must be importable; some scripts use `sys.path.append`)

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
- **Color scheme** ("Haruhi theme" in pysam variants):
  - UTR: `#f50002` (red)
  - CDS: `#53a1b5` (teal)
  - Forward strand reads: `#5d3c2d` (brown)
  - Reverse strand reads: `#f5b024` (gold)
  - Original drawer uses: CDS=yellow, UTR=blue, reads=pink

## Key Conventions for AI Assistants

- **Python 2 codebase**: Use `from __future__ import print_function`. Avoid Python 3-only syntax. The `.sort()` calls use the old pandas API (not `.sort_values()`).
- **No test suite or CI/CD**: There are no automated tests. Verify changes manually.
- **No package manager**: Dependencies are installed globally via pip2. No requirements.txt or setup.py exists.
- **Hardcoded paths**: Many pysam-based scripts contain hardcoded absolute paths to reference files and BAM files (e.g., `/ref/analysis/...`, `/mnt/c/ubuntu.download/...`). These must be edited per-environment.
- **Script naming**: Files are numbered by pipeline step (`1.` = indexing, `2.` = drawing). Variants are appended as dot-separated descriptors.
- **Large binary files in repo**: GFF3 files, pickle files, SAM files, and images are committed directly. Be mindful of repo size.
- **Coordinate system**: GFF3 uses 1-based coordinates. A 100bp padding is added on each side of the gene region for visualization.
- **CIGAR parsing**: Uses regex `r'(\d+)(\w)'` to parse CIGAR strings. Handles M (match), N (intron skip), and I (insertion) operations.
- **Read filtering**: Skips improperly paired reads, duplicates, QC-failed reads, and secondary alignments.
