#!/usr/bin/env python
"""
BAMVIS-GENE: RNA-seq BAM alignment visualization at the gene level.

Unified CLI tool that consolidates the two-step pipeline:
  1. bamvis.py index  — Parse GFF3 into a pandas pickle index
  2. bamvis.py draw   — Render SVG from pickle index + BAM file(s)

Usage:
  pip install -r requirements.txt
  python bamvis.py index  <gff3_file> [--no-longest]
  python bamvis.py draw   --pickle <pk> --bam <bam> --gene <gene>
  python bamvis.py draw   --pickle <pk> --bam <bam> --region chr:start-end
"""

from __future__ import print_function
import sys
import os
import argparse
import re
import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Utility: kang.py functions inlined for portability
# ---------------------------------------------------------------------------

def Fasta2dic(file_fasta):
    dic = {}
    bulk = open(file_fasta).read()
    for each in bulk.split('>'):
        if each == '':
            continue
        strHD = each.split('\n')[0].split()[0]
        strSeq = ''.join(each.split('\n')[1:])
        dic[strHD] = strSeq
    return dic

def rev_comp(strSeq):
    dicComp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    strCseq = ''
    for i in strSeq:
        try:
            strCseq += dicComp[i.upper()]
        except KeyError:
            strCseq += 'N'
    return strCseq[::-1]

# ---------------------------------------------------------------------------
# Step 1: GFF3 Indexing
# ---------------------------------------------------------------------------

def cmd_index(args):
    """Parse a GFF3 file into a pandas pickle index."""
    file_gff = args.gff3
    if not os.path.isfile(file_gff):
        print("Error: GFF3 file not found: %s" % file_gff, file=sys.stderr)
        sys.exit(1)

    print("Reading GFF3: %s" % file_gff)
    df_gff = pd.read_csv(file_gff, sep='\t', comment='#', header=None)

    df_gff_mRNA = df_gff[(df_gff[2] == 'mRNA')].copy()

    def infoparse(x):
        key = [i.split('=')[0] for i in x.split(';')]
        value = [i.split('=')[1] for i in x.split(';')]
        return dict(zip(key, value))

    def get_genename(x):
        if x[2] == 'gene':
            return None
        dic = infoparse(x[8])
        return dic['Parent']

    def get_transcript(x):
        dic = infoparse(x[8])
        if x[2] == 'gene':
            return None
        elif x[2] == 'mRNA':
            return dic['ID']
        else:
            return dic['Parent']

    df_gff_mRNA['genename'] = df_gff_mRNA.apply(get_genename, axis=1)
    df_gff_mRNA['transcript'] = df_gff_mRNA.apply(get_transcript, axis=1)
    df_gff_mRNA_ix = df_gff_mRNA.set_index('transcript')
    df_gff['transcriptname'] = df_gff.apply(get_transcript, axis=1)

    def get_genename2(x):
        if x is not None:
            try:
                return df_gff_mRNA_ix.loc[x]['genename']
            except KeyError:
                return None
        else:
            return None

    df_gff['genename'] = df_gff['transcriptname'].apply(get_genename2)

    outfile = file_gff + '.pandas.df.pk'

    if args.no_longest:
        # Index by genename only, keep all transcripts
        df_gff_index = df_gff.set_index(['genename'])
        print("Indexing by genename (all transcripts)")
    else:
        # Index by genename + longest flag
        try:
            df_gff_mRNA_ix['longest'] = df_gff_mRNA_ix[8].apply(lambda x: infoparse(x)['longest'])
        except KeyError:
            print("Warning: 'longest' attribute not found in GFF3. Falling back to --no-longest mode.",
                  file=sys.stderr)
            df_gff_index = df_gff.set_index(['genename'])
            df_gff_index.to_pickle(outfile)
            print("Saved: %s" % outfile)
            return

        def get_longest(x):
            if x is None:
                return None
            try:
                return df_gff_mRNA_ix.loc[x]['longest']
            except KeyError:
                return None

        df_gff['longest'] = df_gff['transcriptname'].apply(get_longest)
        df_gff_index = df_gff.set_index(['genename', 'longest'])
        print("Indexing by genename + longest transcript")

    df_gff_index.to_pickle(outfile)
    print("Saved: %s" % outfile)


# ---------------------------------------------------------------------------
# Step 2: BAM Visualization (SVG Drawing)
# ---------------------------------------------------------------------------

# SVG templates
RECT = '<rect x="%d" y="%d" width="%d" height="%d" style="fill:%s;stroke:%s;stroke-width:1;fill-opacity:1;stroke-opacity:0.3" />'
SVG_LINE = '<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="%d" />'
TEXT = '<text x="%d"  y="%d" style="font-family:Arial;font-size:%d;stroke:black;fill:black;">%s</text>'

# Color scheme (Haruhi theme)
C_UTR = '#f50002'
C_CDS = '#53a1b5'
C_FSTRAND = '#5d3c2d'
C_SSTRAND = '#f5b024'

# Layout constants
EACH_HEIGHT = 5
EACH_SPACE = 2
GENESTART = 3
CANVAS_WIDTH = 1000


def cigar_parse(cigar):
    return re.findall(r'(\d+)(\w)', cigar)


def get_ratio(x, real_width):
    return int(CANVAS_WIDTH * float(x) / float(real_width))


def init_svg(f, canvas_width, canvas_height):
    print('''<?xml version="1.0" encoding="utf-8" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN"
  "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
<svg height="%s" width="%s" version="1.1" viewBox="0 0 %s %s" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">''' % (
        canvas_height, canvas_width, canvas_width, canvas_height), file=f)


def end_svg(f):
    print('</svg>', file=f)


def draw_words(y, strWord, f):
    print(TEXT % (15, y + 8, 10, strWord), file=f)
    return y + 8 + EACH_SPACE + 5


def draw_gene(start_height, left, right, genename, chromosome, strand, df_gff_ix, index_type, f):
    """Draw gene model (mRNA line + CDS/UTR rectangles)."""
    real_width = right - left + 1 + 200

    if index_type == 'longest':
        df = df_gff_ix.loc[(genename, '1')]
    elif index_type == 'transcript':
        df = df_gff_ix.loc[(genename)]
    else:
        df = df_gff_ix.loc[(genename)]

    df = df.reset_index()
    df = df.sort_values(by=3)

    for i in df.index:
        if df.loc[i][2] == 'mRNA':
            strand = df.loc[i][6]
            mx1 = get_ratio(df.loc[i][3] - left + 100, real_width)
            my1 = start_height + int(EACH_HEIGHT / 2) + GENESTART + 1
            mx2 = get_ratio(df.loc[i][4] - left + 100, real_width)
            my2 = start_height + int(EACH_HEIGHT / 2) + GENESTART + 1
            print(SVG_LINE % (mx1, my1, mx2, my2, 'black', 1), file=f)
    for i in df.index:
        if 'UTR' in str(df.loc[i][2]):
            ux1 = get_ratio(df.loc[i][3] - left + 100, real_width)
            uy1 = start_height + GENESTART + 1
            uw = get_ratio(df.loc[i][4] - df.loc[i][3], real_width)
            uh = EACH_HEIGHT
            print(RECT % (ux1, uy1, uw, uh, C_UTR, 'black'), file=f)
        elif 'CDS' in str(df.loc[i][2]):
            cx1 = get_ratio(df.loc[i][3] - left + 100, real_width)
            cy1 = start_height + GENESTART + 1
            cw = get_ratio(df.loc[i][4] - df.loc[i][3], real_width)
            ch = EACH_HEIGHT
            print(RECT % (cx1, cy1, cw, ch, C_CDS, 'black'), file=f)
    print(TEXT % (mx1, GENESTART + EACH_HEIGHT * 4, 10,
                  '%s (%s) %s:%d-%d' % (genename, strand, chromosome, left, right)), file=f)
    return start_height + GENESTART + EACH_HEIGHT * 4


def draw_alignment(start_height, chromosome, left, right, file_bam, f, total_canvas_rows=20):
    """Draw aligned reads from a BAM file."""
    real_width = right - left + 1 + 200
    gene_space = np.zeros([total_canvas_rows, real_width])
    samfile = pysam.AlignmentFile(file_bam, "rb")
    it = samfile.fetch(chromosome, int(left), int(right))

    dicFID2row = {}
    dicFID2line = {}

    # First pass: collect reads by fragment (query_name)
    for line in tqdm(it, desc="  Reading %s" % os.path.basename(file_bam)):
        if line.is_paired and not line.is_proper_pair:
            continue
        if line.is_duplicate:
            continue
        if line.is_qcfail:
            continue
        if line.is_secondary:
            continue
        try:
            dicFID2line[line.query_name].append(line)
        except KeyError:
            dicFID2line[line.query_name] = [line]

    dicFID2line_list = sorted(dicFID2line.keys(),
                              key=lambda x: int(dicFID2line[x][0].reference_start))

    # Second pass: pack and draw reads
    for fragmentid in tqdm(dicFID2line_list, desc="  Drawing reads"):
        lines = sorted(dicFID2line[fragmentid], key=lambda x: int(x.reference_start))
        line = lines[0]
        fragstart = line.reference_start - left + 1 + 100
        if fragstart < 0:
            continue

        cigar = line.cigarstring
        try:
            cigarM = cigar_parse(cigar)
        except TypeError:
            continue
        cigarstrings = [x[1] for x in cigarM]
        cigarvalues = [x[0] for x in cigarM]

        if line.is_paired:
            if not line.is_proper_pair:
                continue
            if len(lines) < 2:
                continue
            try:
                fragmentsize = line.mpos + lines[1].reference_length - line.reference_start
            except TypeError:
                continue
            bridge_length = line.mpos - line.reference_start - line.reference_length
            if fragmentsize > 0 and bridge_length < 10:
                continue
        else:
            fragmentsize = line.reference_length

        if fragmentsize is None or abs(fragmentsize) == 0:
            continue
        if fragmentsize < 0:
            continue

        # Row selection (greedy packing)
        try:
            srow = dicFID2row[fragmentid]
        except KeyError:
            srow = None
            for nrow, row in enumerate(gene_space):
                end_idx = min(fragstart + fragmentsize, real_width)
                if fragstart >= real_width or end_idx <= fragstart:
                    break
                if max(row[fragstart:end_idx]) > 0:
                    continue
                else:
                    srow = nrow
                    gene_space[srow, fragstart:end_idx] += 1
                    break
            dicFID2row[fragmentid] = srow

        if srow is None:
            continue

        if line.is_paired:
            if fragstart < 0:
                continue
            # Determine strand color
            if (line.is_read1 and not line.is_reverse) or \
               (line.is_read2 and line.is_reverse):
                blockcolor = C_SSTRAND
            else:
                blockcolor = C_FSTRAND

            # Draw V-shaped bridge between mates
            coverleng = 0
            for j, cs in enumerate(cigarstrings):
                if cs in ('M', 'N', 'I'):
                    coverleng += int(cigarvalues[j])
            x1 = get_ratio(fragstart + coverleng, real_width)
            y1 = start_height + srow * (EACH_HEIGHT + EACH_SPACE) + int(EACH_HEIGHT / 2)
            y2 = y1 + int(EACH_HEIGHT / 2)
            x3 = get_ratio(line.mpos - left + 1 + 100, real_width)
            x2 = (x3 + x1) // 2
            print(SVG_LINE % (x1, y1, x2, y2, 'black', 1), file=f)
            print(SVG_LINE % (x2, y2, x3, y1, 'black', 1), file=f)

            # Draw each mate's CIGAR blocks
            for eline in lines:
                w = 0
                efragstart = eline.reference_start - left + 1 + 100
                ecigar = eline.cigarstring
                try:
                    ecigarM = cigar_parse(ecigar)
                except TypeError:
                    continue
                ecigarstrings = [x[1] for x in ecigarM]
                ecigarvalues = [x[0] for x in ecigarM]

                for j, cs in enumerate(ecigarstrings):
                    if cs == 'M':
                        x = get_ratio(w + efragstart, real_width)
                        if x < 0:
                            continue
                        y = start_height + srow * (EACH_HEIGHT + EACH_SPACE)
                        ew = int(ecigarvalues[j])
                        w += ew
                        print(RECT % (x, y, get_ratio(ew, real_width), EACH_HEIGHT, blockcolor, 'black'), file=f)
                    elif cs in ('N', 'I'):
                        x1 = get_ratio(w + efragstart, real_width)
                        if x1 < 0:
                            continue
                        y1 = start_height + srow * (EACH_HEIGHT + EACH_SPACE) + int(EACH_HEIGHT / 2)
                        ew = int(ecigarvalues[j])
                        w += ew
                        x2 = get_ratio(w + efragstart, real_width)
                        print(SVG_LINE % (x1, y1, x2, y1, 'pink', 1), file=f)

        else:
            # Single-end read
            fragstart = line.reference_start - left + 1 + 100
            if fragstart < 0:
                continue
            blockcolor = C_SSTRAND if line.is_reverse else C_FSTRAND
            w = 0
            for j, cs in enumerate(cigarstrings):
                if cs == 'M':
                    x = get_ratio(w + fragstart, real_width)
                    if x < 0:
                        continue
                    y = start_height + srow * (EACH_HEIGHT + EACH_SPACE)
                    ew = int(cigarvalues[j])
                    w += ew
                    print(RECT % (x, y, get_ratio(ew, real_width), EACH_HEIGHT, blockcolor, 'black'), file=f)
                elif cs in ('N', 'I'):
                    x1 = get_ratio(w + fragstart, real_width)
                    if x1 < 0:
                        continue
                    y1 = start_height + srow * (EACH_HEIGHT + EACH_SPACE) + int(EACH_HEIGHT / 2)
                    ew = int(cigarvalues[j])
                    w += ew
                    x2 = get_ratio(w + fragstart, real_width)
                    print(SVG_LINE % (x1, y1, x2, y1, 'pink', 1), file=f)

    samfile.close()
    return start_height + (EACH_HEIGHT + EACH_SPACE) * len(gene_space)


def cmd_draw(args):
    """Generate SVG visualization from pickle index + BAM file(s)."""
    pickle_file = args.pickle
    bam_files = args.bam
    total_rows = args.rows
    output = args.output

    if not os.path.isfile(pickle_file):
        print("Error: Pickle file not found: %s" % pickle_file, file=sys.stderr)
        sys.exit(1)
    for bam in bam_files:
        if not os.path.isfile(bam):
            print("Error: BAM file not found: %s" % bam, file=sys.stderr)
            sys.exit(1)

    print("Loading pickle index: %s" % pickle_file)
    df_gff_index = pd.read_pickle(pickle_file)

    # Detect index type from pickle columns
    index_names = list(df_gff_index.index.names)
    if 'longest' in index_names:
        index_type = 'longest'
        df_gff_ix = df_gff_index.reset_index().set_index(['genename', 'longest', 2])
    elif 'transcriptname' in df_gff_index.columns:
        index_type = 'transcript'
        df_gff_ix = df_gff_index.reset_index().set_index(['transcriptname', 2])
    else:
        index_type = 'genename'
        df_gff_ix = df_gff_index.reset_index().set_index(['genename', 2])
    df_gff_ix.sort_index(inplace=True)

    # --- Gene mode ---
    if args.gene:
        genename = args.gene
        print("Gene mode: %s" % genename)

        if index_type == 'longest':
            try:
                strand = df_gff_ix.loc[(genename, '1', 'mRNA')][6].values[0]
                chromosome, left, right = df_gff_ix.loc[(genename, '1', 'mRNA')][[0, 3, 4]].values[0]
            except KeyError:
                print("Error: Gene '%s' not found in index." % genename, file=sys.stderr)
                sys.exit(1)
        else:
            try:
                strand = df_gff_ix.loc[(genename, 'mRNA')][6].values[0]
                chromosome, left, right = df_gff_ix.loc[(genename, 'mRNA')][[0, 3, 4]].values[0]
            except KeyError:
                print("Error: Gene '%s' not found in index." % genename, file=sys.stderr)
                sys.exit(1)

        left = int(left)
        right = int(right)

        if output is None:
            output = '%s.bamvisgene.svg' % genename

        real_width = right - left + 1 + 200
        canvas_height = (total_rows + 10) * (EACH_HEIGHT + EACH_SPACE) * len(bam_files) + 1000

        print("Region: %s:%d-%d (%s)" % (chromosome, left, right, strand))
        print("Output: %s" % output)

        with open(output, 'w') as f:
            init_svg(f, CANVAS_WIDTH, canvas_height)
            endp = draw_gene(5, left, right, genename, chromosome, strand, df_gff_ix, index_type, f)
            for bam in bam_files:
                endp = draw_alignment(endp, chromosome, left, right, bam, f, total_rows)
                endp = draw_words(endp, os.path.basename(bam), f)
            end_svg(f)

        print("Done: %s" % output)

    # --- Region mode ---
    elif args.region:
        region_str = args.region
        # Parse chr:start-end
        m = re.match(r'^(.+):(\d+)-(\d+)$', region_str)
        if not m:
            print("Error: Invalid region format. Use chr:start-end (e.g. Chr1:1000-5000)", file=sys.stderr)
            sys.exit(1)
        chromosome = m.group(1)
        left = int(m.group(2))
        right = int(m.group(3))
        strand = args.strand if args.strand else '+'

        genename = '%s_%d_%d' % (chromosome, left, right)

        if output is None:
            output = '%s_%d_%d.bamvisgene.svg' % (chromosome, left, right)

        real_width = right - left + 1 + 200
        canvas_height = (total_rows + 10) * (EACH_HEIGHT + EACH_SPACE) * len(bam_files) + 1000

        print("Region mode: %s:%d-%d (%s)" % (chromosome, left, right, strand))
        print("Output: %s" % output)

        with open(output, 'w') as f:
            init_svg(f, CANVAS_WIDTH, canvas_height)
            endp = 5
            for bam in bam_files:
                endp = draw_alignment(endp, chromosome, left, right, bam, f, total_rows)
                endp = draw_words(endp, os.path.basename(bam), f)
            end_svg(f)

        print("Done: %s" % output)

    else:
        print("Error: Specify --gene <name> or --region chr:start-end", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog='bamvis.py',
        description='BAMVIS-GENE: Visualize RNA-seq BAM alignments as SVG at the gene level.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''Examples:
  # Step 1: Index a GFF3 file
  python bamvis.py index Athaliana_167_TAIR10.gene.gff3
  python bamvis.py index Creinhardtii_281_v5.5.gene.gff3 --no-longest

  # Step 2: Draw a gene
  python bamvis.py draw --pickle Athaliana_167_TAIR10.gene.gff3.pandas.df.pk \\
                        --bam sample1.sorted.bam --bam sample2.sorted.bam \\
                        --gene AT1G01010

  # Step 2: Draw a genomic region (no gene model)
  python bamvis.py draw --pickle Creinhardtii_281_v5.5.gene.gff3.pandas.df.pk \\
                        --bam reads.sorted.bam \\
                        --region Chr9:3586521-3608027 --strand +
''')

    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # --- index subcommand ---
    p_index = subparsers.add_parser('index', help='Parse GFF3 into a pandas pickle index')
    p_index.add_argument('gff3', help='GFF3 annotation file')
    p_index.add_argument('--no-longest', action='store_true',
                         help='Keep all transcripts (do not filter for longest)')

    # --- draw subcommand ---
    p_draw = subparsers.add_parser('draw', help='Render SVG from pickle index + BAM file(s)')
    p_draw.add_argument('--pickle', required=True, help='Pandas pickle index file (.pandas.df.pk)')
    p_draw.add_argument('--bam', action='append', required=True,
                        help='BAM file (repeat --bam for multiple files)')
    p_draw.add_argument('--gene', help='Gene name to visualize')
    p_draw.add_argument('--region', help='Genomic region as chr:start-end (e.g. Chr1:1000-5000)')
    p_draw.add_argument('--strand', choices=['+', '-'], default='+',
                        help='Strand for region mode (default: +)')
    p_draw.add_argument('--rows', type=int, default=20,
                        help='Max read packing rows (default: 20)')
    p_draw.add_argument('--output', '-o', help='Output SVG filename (auto-generated if omitted)')

    args = parser.parse_args()

    if args.command == 'index':
        cmd_index(args)
    elif args.command == 'draw':
        cmd_draw(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == '__main__':
    main()
