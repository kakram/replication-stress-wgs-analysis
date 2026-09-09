#!/usr/bin/env python3
"""
rt_celline_matched.py -- replication-timing analysis with each cell line
measured against its OWN timing track, at matched resolution.

Supersedes the replication-timing section of regional_analysis.py, which
used the MCF-7 track for both cell lines and therefore attenuated any real
U2OS signal.

Tracks
------
U2OS   GEO GSE211592, Repli-seq log2(Early/Late), IDH2-wild-type U2OS,
       three biological replicates, hg19, 50 kb bins. Consensus = mean of
       the three replicates over bins present in all three.
MCF-7  the existing Repli-seq bedGraph, hg19, 1 kb bins, aggregated here
       to the same 50 kb grid so the two cell lines are compared at
       identical resolution.

Both tracks use "higher = earlier replication" but on different scales
(MCF-7 roughly 0-86; U2OS log2 ratio roughly -6 to +6). Raw values are
therefore NOT comparable between cell lines. Every cross-line statement
uses within-track percentile rank instead, which is scale-free.

Variants are lifted GRCh38 -> hg19 to match the tracks, restricted to
primary contigs.
"""
import glob, gzip, json, os, re, sys

import numpy as np
from pyliftover import LiftOver
from scipy.stats import mannwhitneyu, spearmanr, pearsonr

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'results', 'regional')
RT_DIR = os.path.join(DATA, 'annotation', 'U2OS_RT')
MCF7_RT = os.path.join(DATA, 'MCF-7-RT-BedGraph.bedgrpah')
CHAIN = os.path.join(DATA, 'annotation', 'hg38ToHg19.over.chain.gz')

BIN = 50000
MIN_DP, MIN_VAF, MIN_ALT, MIN_POPAF = 20, 0.05, 3, 3.0
PRIMARY = {f'chr{i}' for i in range(1, 23)} | {'chrX'}
POPAF_RE = re.compile(r'(?:^|;)POPAF=([^;]+)')
RNG = np.random.default_rng(20260904)

CELL_LINE = {'WTUN': 'MCF7', 'WTAPH': 'MCF7', 'B4UN': 'MCF7', 'B4APH': 'MCF7',
             'WT-U_cleaned': 'U2OS', 'WT-A_cleaned': 'U2OS',
             'G3-U_cleaned': 'U2OS', 'G3-A_cleaned': 'U2OS'}
ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']


def read_bedgraph(path, opener=open):
    """-> dict[(chrom, bin_index)] = value, at BIN resolution."""
    d = {}
    with opener(path, 'rt') as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 4 or p[0] not in PRIMARY:
                continue
            d[(p[0], int(p[1]) // BIN)] = float(p[3])
    return d


def build_u2os():
    reps = sorted(glob.glob(os.path.join(RT_DIR, '*_IDH2_W_rep*.bedgraph.gz')))
    if len(reps) < 2:
        sys.exit(f'Expected >=2 wild-type U2OS replicates, found {len(reps)}')
    tables = [read_bedgraph(p, gzip.open) for p in reps]
    shared = set(tables[0])
    for t in tables[1:]:
        shared &= set(t)
    shared = sorted(shared)
    mat = np.array([[t[k] for k in shared] for t in tables])

    print(f'U2OS: {len(reps)} wild-type replicates, {len(shared):,} shared 50 kb bins')
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            r, _ = pearsonr(mat[i], mat[j])
            rho, _ = spearmanr(mat[i], mat[j])
            print(f'   rep{i+1} vs rep{j+1}:  Pearson r={r:.3f}  Spearman rho={rho:.3f}')
    return {k: v for k, v in zip(shared, mat.mean(axis=0))}, len(reps), \
           [[float(pearsonr(mat[i], mat[j])[0]) for j in range(len(reps))]
            for i in range(len(reps))]


def build_mcf7():
    """Aggregate the 1 kb MCF-7 track to the 50 kb grid."""
    acc = {}
    with open(MCF7_RT) as fh:
        for line in fh:
            p = line.split('\t')
            if len(p) < 4 or p[0] not in PRIMARY:
                continue
            k = (p[0], int(p[1]) // BIN)
            s, n = acc.get(k, (0.0, 0))
            acc[k] = (s + float(p[3]), n + 1)
    print(f'MCF-7: 1 kb track aggregated to {len(acc):,} bins of 50 kb')
    return {k: s / n for k, (s, n) in acc.items()}


def percentiles(track):
    """value -> percentile rank within the track (0-100)."""
    vals = np.array(sorted(track.values()))
    return vals


def variants(vcf_path):
    paired = None
    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                paired = len(line.rstrip('\n').split('\t')) > 10
                continue
            c = line.rstrip('\n').split('\t')
            if c[6] != 'PASS' or c[0] not in PRIMARY:
                continue
            keys = c[8].split(':')
            t = dict(zip(keys, c[9].split(':')))
            try:
                if int(t.get('DP', 0)) < MIN_DP: continue
                if int(t['AD'].split(',')[1]) < MIN_ALT: continue
                if float(t['AF'].split(',')[0]) < MIN_VAF: continue
            except (KeyError, ValueError, IndexError):
                continue
            if paired:
                n = dict(zip(keys, c[10].split(':')))
                try:
                    if int(n.get('DP', 0)) < MIN_DP: continue
                    if int(n['AD'].split(',')[1]) != 0: continue
                except (KeyError, ValueError, IndexError):
                    continue
            m = POPAF_RE.search(c[7])
            if m and float(m.group(1).split(',')[0]) < MIN_POPAF:
                continue
            yield c[0], int(c[1])


def main():
    os.makedirs(OUT, exist_ok=True)
    u2os, n_reps, corr = build_u2os()
    mcf7 = build_mcf7()
    tracks = {'U2OS': u2os, 'MCF7': mcf7}
    sorted_vals = {k: percentiles(v) for k, v in tracks.items()}
    lo = LiftOver(CHAIN)

    # background: all bins of each track, as percentiles (uniform by construction)
    bg = {k: np.array(sorted(v.values())) for k, v in tracks.items()}

    res = {'u2os_replicates': n_reps, 'u2os_replicate_pearson': corr,
           'bin_bp': BIN, 'samples': {}}

    print()
    for path in sorted(glob.glob(os.path.join(DATA, '*', '*', '*_somatic.vcf.gz'))):
        name = os.path.basename(path).replace('_somatic.vcf.gz', '')
        line = CELL_LINE.get(name)
        if line is None:
            continue
        track, svals = tracks[line], sorted_vals[line]
        vals = []
        for chrom, pos in variants(path):
            conv = lo.convert_coordinate(chrom, pos - 1)
            if not conv or conv[0][0] != chrom:
                continue
            v = track.get((chrom, conv[0][1] // BIN))
            if v is not None:
                vals.append(v)
        vals = np.array(vals)
        pct = 100 * np.searchsorted(svals, vals) / len(svals)
        u, p = mannwhitneyu(vals, bg[line], alternative='two-sided')
        rbc = 1 - 2 * u / (len(vals) * len(bg[line]))
        res['samples'][name] = {
            'cell_line': line, 'n': int(len(vals)),
            'rt_mean': round(float(vals.mean()), 4),
            'rt_percentile_mean': round(float(pct.mean()), 2),
            'rt_percentile_median': round(float(np.median(pct)), 2),
            'effect_vs_background': round(float(rbc), 4),
            'p': float(p),
            'pct_in_latest_quintile': round(float((pct < 20).mean() * 100), 2),
        }
        print(f'  {name:<15} {line:<5} n={len(vals):>7,} '
              f'RTpct={pct.mean():>6.2f}  effect={rbc:>7.4f}  p={p:.2e}')

    with open(os.path.join(OUT, 'rt_celline_matched.json'), 'w') as fh:
        json.dump(res, fh, indent=2)

    print('\n=== each cell line against its OWN track, 50 kb bins ===')
    print(f"{'sample':<16}{'line':<7}{'n':>9}{'RT %ile':>10}"
          f"{'effect':>10}{'% in latest 20%':>18}")
    print('-' * 70)
    for s in ORDER:
        r = res['samples'].get(s)
        if not r: continue
        print(f"{s:<16}{r['cell_line']:<7}{r['n']:>9,}"
              f"{r['rt_percentile_mean']:>10}{r['effect_vs_background']:>10}"
              f"{r['pct_in_latest_quintile']:>17}%")
    print('\n(background mean percentile is 50.0 by construction; lower = later '
          'replicating. "% in latest 20%" against a 20% null.)')


if __name__ == '__main__':
    main()
