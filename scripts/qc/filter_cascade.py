#!/usr/bin/env python3
"""
filter_cascade.py -- tiered filtering of the Azenta somatic call sets.

Applies a cumulative filter cascade to each call set and reports variant
counts, composition and Ts/Tv at every tier. The point is to establish how
much of each call set survives progressively stricter criteria, and whether
Ts/Tv rises toward the ~2.0 expected of genuine variants as it does.

Tiers
-----
  T0  ALL     every record in the VCF
  T1  PASS    caller FILTER == PASS
  T2  DEPTH   T1 + DP >= 20 in the tumour AND, for paired call sets, DP >= 20
              in the matched normal. This equalises calling power across
              contrasts that were sequenced to different depths.
  T3  VAF     T2 + tumour allele fraction >= 0.05 and >= 3 supporting reads
  T4  POP     T3 + POPAF >= 3, i.e. population allele frequency <= 1e-3
              (Mutect2/TNhaplotyper2 phred-scaled gnomAD frequency)
  T5  PRIVATE T4 + zero alt-supporting reads in the matched normal
              (paired call sets only; strict knockout-private set)

Also records, for paired call sets, the distribution of alt-read support in
the matched normal among T4 variants -- the clonal-bottleneck diagnostic.
A variant with alt reads in the wild-type at good depth pre-existed the
knockout and cannot be attributed to loss of the gene.

Usage
-----
  python3 scripts/qc/filter_cascade.py <sample.vcf.gz>
  python3 scripts/qc/filter_cascade.py --all      # any not yet done
  python3 scripts/qc/filter_cascade.py --table
"""
import argparse, glob, gzip, json, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'results', 'cascade')

TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
TIERS = ['T0_ALL', 'T1_PASS', 'T2_DEPTH', 'T3_VAF', 'T4_POP', 'T5_PRIVATE']

MIN_DP = 20
MIN_VAF = 0.05
MIN_ALT_READS = 3
MIN_POPAF = 3.0          # population AF <= 1e-3


def new_bucket():
    return dict(n=0, snv=0, ins=0, dele=0, mnv=0, ts=0, tv=0)


def classify(bucket, ref, alt):
    bucket['n'] += 1
    if len(ref) == 1 and len(alt) == 1:
        bucket['snv'] += 1
        if (ref, alt) in TRANSITIONS:
            bucket['ts'] += 1
        else:
            bucket['tv'] += 1
    elif len(ref) == len(alt):
        bucket['mnv'] += 1
    elif len(alt) > len(ref):
        bucket['ins'] += 1
    else:
        bucket['dele'] += 1


def parse_sample(fmt_keys, field):
    """Return (dp, alt_reads, vaf) for one sample column."""
    vals = field.split(':')
    d = dict(zip(fmt_keys, vals))
    try:
        dp = int(d.get('DP', '0'))
    except ValueError:
        dp = 0
    alt = 0
    ad = d.get('AD')
    if ad:
        parts = ad.split(',')
        if len(parts) > 1:
            try:
                alt = int(parts[1])
            except ValueError:
                alt = 0
    try:
        vaf = float(d.get('AF', '0').split(',')[0])
    except ValueError:
        vaf = 0.0
    return dp, alt, vaf


POPAF_RE = re.compile(r'(?:^|;)POPAF=([^;]+)')


def run(vcf_path):
    name = os.path.basename(vcf_path).replace('_somatic.vcf.gz', '')
    dest = os.path.join(OUT, f'{name}_cascade.json')

    buckets = {t: new_bucket() for t in TIERS}
    normal_support = {'0': 0, '1': 0, '2': 0, '3+': 0}
    paired = None
    dropped_by_depth = 0

    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                paired = len(line.rstrip('\n').split('\t')) > 10
                continue
            c = line.rstrip('\n').split('\t')
            ref, alt, filt, info, fmt = c[3], c[4].split(',')[0], c[6], c[7], c[8]
            fmt_keys = fmt.split(':')

            classify(buckets['T0_ALL'], ref, alt)
            if filt != 'PASS':
                continue
            classify(buckets['T1_PASS'], ref, alt)

            t_dp, t_alt, t_vaf = parse_sample(fmt_keys, c[9])
            n_dp = n_alt = 0
            if paired:
                n_dp, n_alt, _ = parse_sample(fmt_keys, c[10])

            if t_dp < MIN_DP or (paired and n_dp < MIN_DP):
                dropped_by_depth += 1
                continue
            classify(buckets['T2_DEPTH'], ref, alt)

            if t_vaf < MIN_VAF or t_alt < MIN_ALT_READS:
                continue
            classify(buckets['T3_VAF'], ref, alt)

            m = POPAF_RE.search(info)
            popaf = float(m.group(1).split(',')[0]) if m else 99.0
            if popaf < MIN_POPAF:
                continue
            classify(buckets['T4_POP'], ref, alt)

            if paired:
                key = '0' if n_alt == 0 else ('1' if n_alt == 1
                                              else ('2' if n_alt == 2 else '3+'))
                normal_support[key] += 1
                if n_alt == 0:
                    classify(buckets['T5_PRIVATE'], ref, alt)

    for b in buckets.values():
        b['tstv'] = round(b['ts'] / b['tv'], 3) if b['tv'] else None
        b['indel_frac'] = (round((b['ins'] + b['dele']) / b['n'], 4)
                           if b['n'] else None)

    rec = {
        'sample': name,
        'paired': paired,
        'tiers': buckets,
        'pass_dropped_by_depth': dropped_by_depth,
        'normal_alt_support_at_T4': normal_support if paired else None,
        'thresholds': dict(min_dp=MIN_DP, min_vaf=MIN_VAF,
                           min_alt_reads=MIN_ALT_READS, min_popaf=MIN_POPAF),
    }
    os.makedirs(OUT, exist_ok=True)
    with open(dest, 'w') as fh:
        json.dump(rec, fh, indent=2)
    t = buckets
    print(f"  [done] {name:<14} ALL={t['T0_ALL']['n']:>9,}  PASS={t['T1_PASS']['n']:>8,}"
          f"  DEPTH={t['T2_DEPTH']['n']:>8,}  VAF={t['T3_VAF']['n']:>8,}"
          f"  POP={t['T4_POP']['n']:>8,}  PRIV={t['T5_PRIVATE']['n']:>8,}", flush=True)
    return rec


ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']


def table():
    recs = []
    for p in glob.glob(os.path.join(OUT, '*_cascade.json')):
        with open(p) as fh:
            recs.append(json.load(fh))
    if not recs:
        print('Nothing computed yet.')
        return
    recs.sort(key=lambda r: ORDER.index(r['sample'])
              if r['sample'] in ORDER else 99)

    print('\n=== VARIANT COUNT BY TIER ===')
    hdr = f"{'sample':<15}" + ''.join(f'{t.split("_")[1]:>11}' for t in TIERS)
    print(hdr); print('-' * len(hdr))
    for r in recs:
        print(f"{r['sample']:<15}" +
              ''.join(f"{r['tiers'][t]['n']:>11,}" for t in TIERS))

    print('\n=== Ts/Tv BY TIER ===')
    print(hdr); print('-' * len(hdr))
    for r in recs:
        row = ''
        for t in TIERS:
            v = r['tiers'][t]['tstv']
            row += f"{v if v is not None else '-':>11}"
        print(f"{r['sample']:<15}" + row)

    print('\n=== INDEL FRACTION BY TIER ===')
    print(hdr); print('-' * len(hdr))
    for r in recs:
        row = ''
        for t in TIERS:
            v = r['tiers'][t]['indel_frac']
            row += f"{v if v is not None else '-':>11}"
        print(f"{r['sample']:<15}" + row)

    print('\n=== CLONAL BOTTLENECK: alt reads in matched normal at T4 ===')
    print(f"{'sample':<15}{'0':>10}{'1':>8}{'2':>8}{'3+':>8}{'% with support':>16}")
    for r in recs:
        ns = r.get('normal_alt_support_at_T4')
        if not ns:
            continue
        tot = sum(ns.values())
        sup = tot - ns['0']
        print(f"{r['sample']:<15}{ns['0']:>10,}{ns['1']:>8,}{ns['2']:>8,}"
              f"{ns['3+']:>8,}{100*sup/tot if tot else 0:>15.2f}%")


def pending():
    out = []
    for p in sorted(glob.glob(os.path.join(DATA, '*', '*', '*_somatic.vcf.gz'))):
        name = os.path.basename(p).replace('_somatic.vcf.gz', '')
        if not os.path.exists(os.path.join(OUT, f'{name}_cascade.json')):
            out.append(p)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('vcf', nargs='?')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--table', action='store_true')
    ap.add_argument('--limit', type=int, default=99)
    a = ap.parse_args()
    if a.table:
        table()
    elif a.all:
        for p in pending()[:a.limit]:
            run(p)
    elif a.vcf:
        run(a.vcf)
