#!/usr/bin/env python3
"""
bam_metrics.py -- alignment and coverage metrics for the Azenta BAMs.

Design
------
Read totals come from the BAM index (samtools idxstats) and are EXACT.
Coverage, duplicate rate and proper-pair rate are estimated from a uniform
random sample of windows across chr1-22 and chrX, weighted by chromosome
length. Sampling avoids a full pass over a ~30 GB BAM, which is impractical
over USB, while still giving very tight estimates: a 10 Mb sample at ~25x
contains several million reads, so the sampling error on a rate is well
under 0.1 percentage points (95% CI reported alongside).

Sampling includes assembly gaps, so the 0x figure is comparable to a
whole-genome "% genome no coverage" statistic.

Duplicates, secondary and supplementary alignments and QC-failed reads are
excluded from the coverage calculation (pysam read_callback='all'), so
coverage is post-deduplication.

Usage
-----
  python3 scripts/qc/bam_metrics.py <file.bam>   # one BAM (~1-2 min)
  python3 scripts/qc/bam_metrics.py --next       # next unprocessed in manifest
  python3 scripts/qc/bam_metrics.py --table      # combined table + TSV
"""
import argparse, glob, json, math, os, random, re, statistics, sys, time

import pysam

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, 'results', 'qc')
MANIFEST = os.path.join(REPO, 'scripts', 'qc', 'bam_manifest.txt')
MAIN_CONTIGS = [f'chr{i}' for i in range(1, 23)] + ['chrX']

N_WINDOWS = 2000
WINDOW_BP = 5000
SEED = 20260902


def sample_name(bam_path):
    return re.sub(r'\.aln\.bam$|\.bam$', '', os.path.basename(bam_path))


def wilson_halfwidth(p, n, z=1.96):
    """Half-width of the Wilson 95% interval, in percentage points."""
    if n == 0:
        return None
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    spread = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return round(100 * spread, 4)


def exact_read_totals(bam):
    """Exact mapped/unmapped counts, read from the index only."""
    mapped = unmapped = 0
    for row in bam.get_index_statistics():
        mapped += row.mapped
        unmapped += row.unmapped
    unmapped += bam.nocoordinate
    return mapped, unmapped


def sample_windows(bam, n_windows=N_WINDOWS, window=WINDOW_BP, seed=SEED):
    rng = random.Random(seed)
    lengths = {c: l for c, l in zip(bam.references, bam.lengths)
               if c in MAIN_CONTIGS}
    if not lengths:
        raise SystemExit('No main contigs found')
    contigs = list(lengths)
    weights = [lengths[c] for c in contigs]

    depths = []
    flags = dict(reads=0, dup=0, secondary=0, supplementary=0,
                 paired=0, proper=0, mapq0=0, mapq_sum=0)

    for _ in range(n_windows):
        contig = rng.choices(contigs, weights=weights, k=1)[0]
        limit = lengths[contig] - window
        if limit <= 0:
            continue
        start = rng.randrange(0, limit)
        stop = start + window

        for r in bam.fetch(contig, start, stop):
            if r.reference_start < start:
                continue          # count each read once, in the window it starts in
            flags['reads'] += 1
            if r.is_duplicate:
                flags['dup'] += 1
            if r.is_secondary:
                flags['secondary'] += 1
            if r.is_supplementary:
                flags['supplementary'] += 1
            if r.is_paired:
                flags['paired'] += 1
                if r.is_proper_pair:
                    flags['proper'] += 1
            if r.mapping_quality == 0:
                flags['mapq0'] += 1
            flags['mapq_sum'] += r.mapping_quality

        cols = bam.count_coverage(contig, start, stop,
                                  quality_threshold=0, read_callback='all')
        depths.extend(a + c + g + t for a, c, g, t in zip(*cols))

    return depths, flags


def process(bam_path, force=False):
    name = sample_name(bam_path)
    dest = os.path.join(OUT, f'{name}_bam_metrics.json')
    if os.path.exists(dest) and not force:
        print(f'  [have] {name}', flush=True)
        return None

    t0 = time.time()
    print(f'  [start] {name}', flush=True)
    os.makedirs(OUT, exist_ok=True)

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        mapped, unmapped = exact_read_totals(bam)
        depths, f = sample_windows(bam)

    total = mapped + unmapped
    n = len(depths)
    frac = lambda k: sum(1 for d in depths if d >= k) / n
    rate = lambda a, b: round(100.0 * a / b, 3) if b else None

    rec = {
        'sample': name,
        'bam': bam_path,
        'bam_bytes': os.path.getsize(bam_path),
        # exact, from the index
        'total_reads': total,
        'mapped_reads': mapped,
        'pct_aligned': rate(mapped, total),
        # estimated, from the window sample
        'pct_duplicates': rate(f['dup'], f['reads']),
        'pct_duplicates_ci95': wilson_halfwidth(f['dup'] / f['reads'], f['reads']),
        'pct_properly_paired': rate(f['proper'], f['paired']),
        'pct_properly_paired_ci95': wilson_halfwidth(
            f['proper'] / f['paired'], f['paired']) if f['paired'] else None,
        'pct_mapq0': rate(f['mapq0'], f['reads']),
        'mean_mapq': round(f['mapq_sum'] / f['reads'], 2) if f['reads'] else None,
        'mean_coverage': round(statistics.fmean(depths), 2),
        'median_coverage': statistics.median(depths),
        'pct_genome_0x': round(100 * (1 - frac(1)), 3),
        'pct_genome_10x': round(100 * frac(10), 3),
        'pct_genome_20x': round(100 * frac(20), 3),
        'pct_genome_30x': round(100 * frac(30), 3),
        # provenance
        'reads_sampled': f['reads'],
        'bases_sampled': n,
        'windows': N_WINDOWS,
        'window_bp': WINDOW_BP,
        'seed': SEED,
        'pysam_version': pysam.__version__,
        'runtime_seconds': round(time.time() - t0, 1),
    }
    with open(dest, 'w') as fh:
        json.dump(rec, fh, indent=2)
    print(f'  [done]  {name}  mean={rec["mean_coverage"]}x  '
          f'dup={rec["pct_duplicates"]}%  aligned={rec["pct_aligned"]}%  '
          f'({rec["runtime_seconds"]}s)', flush=True)
    return rec


COLUMNS = [
    ('sample', 'Sample', 16), ('total_reads', 'Total reads', 16),
    ('pct_aligned', '% aligned', 11), ('pct_duplicates', '% dup', 9),
    ('pct_properly_paired', '% proper', 10),
    ('mean_coverage', 'Mean cov', 10), ('median_coverage', 'Median', 8),
    ('pct_genome_0x', '% 0x', 8), ('pct_genome_10x', '% >=10x', 9),
    ('pct_genome_20x', '% >=20x', 9), ('pct_genome_30x', '% >=30x', 9),
]


def table():
    recs = []
    for p in sorted(glob.glob(os.path.join(OUT, '*_bam_metrics.json'))):
        with open(p) as fh:
            recs.append(json.load(fh))
    if not recs:
        print('No metrics yet.')
        return
    order = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
             'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']
    recs.sort(key=lambda r: order.index(r['sample'])
              if r['sample'] in order else 99)
    head = ''.join(f'{lbl:>{w}}' for _, lbl, w in COLUMNS)
    print(head)
    print('-' * len(head))
    for r in recs:
        row = ''
        for k, _, w in COLUMNS:
            v = r.get(k)
            v = f'{v:,}' if isinstance(v, int) and k == 'total_reads' else v
            row += f'{v if v is not None else "NA":>{w}}'
        print(row)
    with open(os.path.join(OUT, 'bam_metrics_summary.tsv'), 'w') as fh:
        fh.write('\t'.join(k for k, _, _ in COLUMNS) + '\n')
        for r in recs:
            fh.write('\t'.join(str(r.get(k, '')) for k, _, _ in COLUMNS) + '\n')
    print(f'\n{len(recs)}/8 samples complete.')


def next_pending():
    with open(MANIFEST) as fh:
        for line in fh:
            bam = line.strip()
            if not bam:
                continue
            if not os.path.exists(os.path.join(
                    OUT, f'{sample_name(bam)}_bam_metrics.json')):
                return bam
    return None


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('bam', nargs='?')
    ap.add_argument('--next', action='store_true')
    ap.add_argument('--table', action='store_true')
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    if a.table:
        table()
    elif a.next:
        b = next_pending()
        if b:
            process(b, a.force)
        else:
            print('All manifest entries processed.')
            table()
    elif a.bam:
        process(a.bam, a.force)
