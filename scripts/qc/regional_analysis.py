#!/usr/bin/env python3
"""
regional_analysis.py -- where do the filtered variants fall?

Tests two positional questions:

  1. Replication timing. Are variants skewed toward late-replicating DNA?
     Uses the MCF-7 Repli-seq track (1 kb bins). Higher signal = EARLIER
     replication on this scale, so a LOWER mean RT value means later.

  2. Common fragile sites. What fraction of variants fall inside annotated
     fragile sites?

Coordinate handling
-------------------
Variants are called on GRCh38; both annotation tracks are hg19. Variant
positions are therefore lifted GRCh38 -> hg19 with the UCSC chain file
(lifting ~400k variants rather than 2.9M annotation bins, and leaving the
reference tracks in their native build). Lift failure rate is reported.

Caveats that must travel with these numbers
-------------------------------------------
* Mutation rate correlates with late replication in essentially all cells.
  A late skew is therefore EXPECTED and is not by itself evidence of any
  ZFP36L1 or aphidicolin effect. The informative comparisons are between
  call sets built by the same design.
* Wild-type sets are tumour-only and knockout sets are paired-subtracted,
  so WT vs KO positional comparisons are design-confounded.
* The fragile-site annotation is the full HumCFS cytogenetic catalogue
  covering ~37% of the genome and mixing rare with aphidicolin-induced
  sites. Absolute "enrichment" against it is uninformative; only relative
  comparisons between samples, which share the annotation, are meaningful.
* Replication timing is from MCF-7. Applying it to U2OS assumes timing is
  conserved between the lines, which is broadly but not perfectly true.
"""
import glob, gzip, json, os, re, sys

import numpy as np
from pyliftover import LiftOver
from scipy.stats import mannwhitneyu

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'results', 'regional')
RT_FILE = os.path.join(DATA, 'MCF-7-RT-BedGraph.bedgrpah')
CFS_DIR = os.path.join(DATA, 'CFS_fixed')
CHAIN = os.path.join(DATA, 'annotation', 'hg38ToHg19.over.chain.gz')

MIN_DP, MIN_VAF, MIN_ALT, MIN_POPAF = 20, 0.05, 3, 3.0
BIN = 1000
POPAF_RE = re.compile(r'(?:^|;)POPAF=([^;]+)')
ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']
RNG = np.random.default_rng(20260904)
PRIMARY = {f'chr{i}' for i in range(1, 23)} | {'chrX'}


def load_rt():
    """chrom -> float array indexed by (pos // 1000); NaN where no data."""
    tracks, maxbin = {}, {}
    with open(RT_FILE) as fh:
        for line in fh:
            c, s, e, v = line.split('\t')
            b = int(s) // BIN
            if b > maxbin.get(c, -1):
                maxbin[c] = b
    for c, m in maxbin.items():
        tracks[c] = np.full(m + 1, np.nan, dtype=np.float32)
    with open(RT_FILE) as fh:
        for line in fh:
            c, s, e, v = line.split('\t')
            tracks[c][int(s) // BIN] = float(v)
    return tracks


def load_cfs(tracks):
    """chrom -> boolean mask over 1 kb bins marking annotated fragile sites.

    Repairs two faults in the delivered BED files: the lowercase `chrx`
    contig name, and FRA10D whose start and end are transposed.
    """
    masks = {c: np.zeros(len(a), dtype=bool) for c, a in tracks.items()}
    n_regions = 0
    for path in sorted(glob.glob(os.path.join(CFS_DIR, '*.bed'))):
        with open(path) as fh:
            for line in fh:
                p = line.rstrip('\n').split('\t')
                if len(p) < 3:
                    continue
                c = 'chrX' if p[0] == 'chrx' else p[0]
                try:
                    s, e = int(p[1]), int(p[2])
                except ValueError:
                    continue
                if s > e:
                    s, e = e, s
                if c not in masks:
                    continue
                lo_b, hi_b = s // BIN, min(e // BIN + 1, len(masks[c]))
                masks[c][lo_b:hi_b] = True
                n_regions += 1
    covered = sum(m.sum() for m in masks.values())
    total = sum(len(m) for m in masks.values())
    print(f'  CFS: {n_regions} regions, {covered * BIN / 1e6:.0f} Mb '
          f'({100 * covered / total:.1f}% of tracked genome)', flush=True)
    return masks


def in_cfs(masks, chrom, pos):
    m = masks.get(chrom)
    if m is None:
        return False
    b = pos // BIN
    return bool(0 <= b < len(m) and m[b])


def variants(vcf_path):
    """Yield hg38 (chrom, pos) for the fully filtered tier."""
    paired = None
    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                paired = len(line.rstrip('\n').split('\t')) > 10
                continue
            c = line.rstrip('\n').split('\t')
            if c[6] != 'PASS':
                continue
            keys = c[8].split(':')
            t = dict(zip(keys, c[9].split(':')))
            try:
                if int(t.get('DP', 0)) < MIN_DP:
                    continue
                if int(t['AD'].split(',')[1]) < MIN_ALT:
                    continue
                if float(t['AF'].split(',')[0]) < MIN_VAF:
                    continue
            except (KeyError, ValueError, IndexError):
                continue
            if paired:
                n = dict(zip(keys, c[10].split(':')))
                try:
                    if int(n.get('DP', 0)) < MIN_DP:
                        continue
                    if int(n['AD'].split(',')[1]) != 0:
                        continue
                except (KeyError, ValueError, IndexError):
                    continue
            m = POPAF_RE.search(c[7])
            if m and float(m.group(1).split(',')[0]) < MIN_POPAF:
                continue
            yield c[0], int(c[1])


def background(tracks, cfs, n=300000):
    """Random positions drawn uniformly over RT-covered genome (hg19)."""
    chroms = [c for c in tracks if c != 'chrY']
    weights = np.array([np.isfinite(tracks[c]).sum() for c in chroms], float)
    weights /= weights.sum()
    picks = RNG.choice(len(chroms), size=n, p=weights)
    rts, incfs = [], 0
    for i in range(n):
        c = chroms[picks[i]]
        b = RNG.integers(0, len(tracks[c]))
        v = tracks[c][b]
        if not np.isfinite(v):
            continue
        rts.append(float(v))
        if in_cfs(cfs, c, b * BIN):
            incfs += 1
    return np.array(rts), incfs / max(len(rts), 1)


def main():
    os.makedirs(OUT, exist_ok=True)
    print('loading tracks...', flush=True)
    tracks = load_rt()
    cfs = load_cfs(tracks)
    lo = LiftOver(CHAIN)

    print('building background...', flush=True)
    bg_rt, bg_cfs_frac = background(tracks, cfs)

    results = {'background': {
        'n': len(bg_rt), 'rt_mean': round(float(bg_rt.mean()), 3),
        'rt_median': round(float(np.median(bg_rt)), 3),
        'cfs_fraction': round(bg_cfs_frac, 4)}}

    for path in sorted(glob.glob(os.path.join(DATA, '*', '*', '*_somatic.vcf.gz'))):
        name = os.path.basename(path).replace('_somatic.vcf.gz', '')
        total = primary = lifted = 0
        rts, cfs_hits, rt_hits = [], 0, 0
        for chrom, pos in variants(path):
            total += 1
            if chrom not in PRIMARY:          # drop unplaced/random/alt/HLA
                continue
            primary += 1
            conv = lo.convert_coordinate(chrom, pos - 1)
            if not conv or conv[0][0] != chrom:
                continue
            lifted += 1
            p19 = conv[0][1]
            if in_cfs(cfs, chrom, p19):
                cfs_hits += 1
            arr = tracks.get(chrom)
            if arr is not None:
                b = p19 // BIN
                if 0 <= b < len(arr) and np.isfinite(arr[b]):
                    rts.append(float(arr[b]))
                    rt_hits += 1
        rts = np.array(rts)
        u, p = mannwhitneyu(rts, bg_rt, alternative='two-sided')
        # rank-biserial correlation = effect size, +ve means later than bg
        rbc = 1 - 2 * u / (len(rts) * len(bg_rt))
        results[name] = {
            'variants_filtered': total,
            'on_primary_contigs': primary,
            'pct_non_primary_contig': round(100 * (1 - primary / total), 3) if total else None,
            'lifted': lifted,
            'pct_primary_unliftable': round(100 * (1 - lifted / primary), 3) if primary else None,
            'rt_n': rt_hits,
            'rt_mean': round(float(rts.mean()), 3),
            'rt_median': round(float(np.median(rts)), 3),
            'rt_vs_background_effect': round(float(rbc), 4),
            'rt_vs_background_p': float(p),
            'cfs_fraction': round(cfs_hits / lifted, 4) if lifted else None,
            'rt_deciles': [round(float(x), 2) for x in
                           np.percentile(rts, np.arange(10, 100, 10))],
        }
        print(f"  {name:<15} n={total:>7,}  nonprim={results[name]['pct_non_primary_contig']}%"
              f"  unliftable={results[name]['pct_primary_unliftable']}%"
              f"  RTmean={results[name]['rt_mean']:>7}"
              f"  CFS={results[name]['cfs_fraction']}", flush=True)

    with open(os.path.join(OUT, 'regional_summary.json'), 'w') as fh:
        json.dump(results, fh, indent=2)

    print('\n=== REPLICATION TIMING (higher = earlier replication) ===')
    print(f"{'set':<16}{'n':>9}{'RT mean':>10}{'RT median':>11}"
          f"{'vs bg':>9}{'p':>12}{'% in CFS':>10}")
    print(f"{'RANDOM BACKGROUND':<16}{results['background']['n']:>9,}"
          f"{results['background']['rt_mean']:>10}"
          f"{results['background']['rt_median']:>11}{'-':>9}{'-':>12}"
          f"{100*results['background']['cfs_fraction']:>9.2f}%")
    print('-' * 77)
    for s in ORDER:
        r = results.get(s)
        if not r:
            continue
        print(f"{s:<16}{r['rt_n']:>9,}{r['rt_mean']:>10}{r['rt_median']:>11}"
              f"{r['rt_vs_background_effect']:>9}{r['rt_vs_background_p']:>12.2e}"
              f"{100*r['cfs_fraction']:>9.2f}%")


if __name__ == '__main__':
    main()
