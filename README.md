# BAMVIS-GENE

RNA-seq BAM alignment visualization at the gene level. Draws gene models (CDS, UTR, mRNA) from GFF3 annotations and overlays aligned sequencing reads from BAM files as SVG images, showing splice junctions, paired-end bridges, and strand orientation.

## Installation

```sh
pip install -r requirements.txt
```

Requires: `pandas`, `pysam`, `numpy`, `tqdm`

## Usage

### Step 1: Index a GFF3 file

```sh
# Filter for longest transcript per gene (default)
python bamvis.py index Athaliana_167_TAIR10.gene.gff3

# Keep all transcripts
python bamvis.py index Creinhardtii_281_v5.5.gene.gff3 --no-longest
```

Output: `<gff3_file>.pandas.df.pk`

### Step 2: Visualize a gene

```sh
python bamvis.py draw --pickle Athaliana_167_TAIR10.gene.gff3.pandas.df.pk \
                      --bam sample.sorted.bam \
                      --gene AT1G01010
```

Multiple BAM files:

```sh
python bamvis.py draw --pickle Athaliana_167_TAIR10.gene.gff3.pandas.df.pk \
                      --bam sample1.sorted.bam \
                      --bam sample2.sorted.bam \
                      --gene AT1G01010
```

### Step 2 (alt): Visualize a genomic region

```sh
python bamvis.py draw --pickle Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk \
                      --bam reads.sorted.bam \
                      --region Chr9:3586521-3608027 --strand +
```

### Options

| Option | Description |
|--------|-------------|
| `--pickle` | Pandas pickle index file (`.pandas.df.pk`) |
| `--bam` | BAM file path (repeat for multiple) |
| `--gene` | Gene name to visualize (draws gene model + reads) |
| `--region` | Genomic region as `chr:start-end` (draws reads only) |
| `--strand` | Strand for region mode (`+` or `-`, default: `+`) |
| `--rows` | Max read packing rows (default: 20) |
| `-o` / `--output` | Output SVG filename (auto-generated if omitted) |

## Output

SVG images (1000px wide) with:
- **Gene model** at the top: mRNA line, CDS rectangles (teal), UTR rectangles (red)
- **Reads** packed vertically: forward strand (brown), reverse strand (gold)
- **Paired-end bridges**: V-shaped lines connecting mate pairs
- **Splice junctions**: thin lines for intron-skipping (CIGAR N) and insertions (CIGAR I)

## Supported organisms

Pre-built pickle indexes are included for:
- *Arabidopsis thaliana* (TAIR10)
- *Chlamydomonas reinhardtii* (v5.5)
- *Eutrema salsugineum* (v1.0)

Any Phytozome-format GFF3 can be indexed with `bamvis.py index`.
