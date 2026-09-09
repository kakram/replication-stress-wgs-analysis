#!/usr/bin/env python3
"""
trinucleotide_context.py -- build the 96-channel mutational matrix.

For every filtered single-base substitution, records the substitution
together with its immediate 5' and 3' neighbours. Substitutions are
normalised to the pyrimidine strand by convention: where the reference
base is A or G, the trinucleotide is reverse-complemented. This gives the
standard 96 categories (6 substitution types x 4 x 4 flanking bases),
which is the input format for COSMIC signature fitting.

Filters match tier T5 of filter_cascade.py for paired call sets and T4 for
tumour-only sets, restricted to primary contigs (chr1-22, X). Centromeric
satellite is not explicitly excluded here, but §D2 of the report shows it
is a substantial component of the tumour-only sets; the `--exclude-unliftable`
option drops positions with no hg19 counterpart, which removes most of it.

Every variant's reference base is checked against the reference genome. A
mismatch rate above a fraction of a percent would indicate an assembly
mismatch, so the rate is reported as a sanity check.

Usage
-----
  python3 scripts/qc/trinucleotide_context.py [--exclude-unliftable]
"""
import argparse, glob, gzip, json, os, re, sys

import pysam

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'results', 'signatures')
FASTA = os.path.expanduser('~/mnt/My Passport for Mac/reference_hg38/hg38.fa')
CHAIN = os.path.join(DATA, 'annotation', 'hg38ToHg19.over.chain.gz')

MIN_DP, MIN_VAF, MIN_ALT, MIN_POPAF = 20, 0.05, 3, 3.0
PRIMARY = {f'chr{i}' for i in range(1, 23)} | {'chrX'}
COMP = str.maketrans('ACGT', 'TGCA')
SUBS = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
BASES = 'ACGT'
CHANNELS = [f'{a}[{s}]{b}' for s in SUBS for a in BASES for b in BASES]
POPAF_RE = re.compile(r'(?:^|;)POPAF=([^;]+)')
ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']


def revcomp(s):
    return s.translate(COMP)[::-1]


def channel(tri, ref, alt):
    """Pyrimidine-normalised 96-channel label, or None if not usable."""
    if len(tri) != 3 or any(b not in BASES for b in tri):
        return None
    if ref in 'AG':
        tri, ref, alt = revcomp(tri), ref.translate(COMP), alt.translate(COMP)
    return f'{tri[0]}[{ref}>{alt}]{tri[2]}'


def filtered_snvs(vcf_path):
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
            ref, alt = c[3], c[4].split(',')[0]
            if len(ref) != 1 or len(alt) != 1 or ref not in BASES or alt not in BASES:
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
            yield c[0], int(c[1]), ref, alt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exclude-unliftable', action='store_true',
                    help='drop positions with no hg19 counterpart (removes '
                         'most GRCh38-only centromeric satellite)')
    args = ap.parse_args()

    fa = pysam.FastaFile(FASTA)
    lo = None
    if args.exclude_unliftable:
        from pyliftover import LiftOver
        lo = LiftOver(CHAIN)

    os.makedirs(OUT, exist_ok=True)
    matrix, meta = {}, {}

    for path in sorted(glob.glob(os.path.join(DATA, '*', '*', '*_somatic.vcf.gz'))):
        name = os.path.basename(path).replace('_somatic.vcf.gz', '')
        counts = {ch: 0 for ch in CHANNELS}
        n = ref_ok = dropped_lift = bad_ctx = 0
        for chrom, pos, ref, alt in filtered_snvs(path):
            n += 1
            if lo is not None:
                conv = lo.convert_coordinate(chrom, pos - 1)
                if not conv or conv[0][0] != chrom:
                    dropped_lift += 1
                    continue
            tri = fa.fetch(chrom, pos - 2, pos + 1).upper()
            if len(tri) == 3 and tri[1] == ref:
                ref_ok += 1
            ch = channel(tri, ref, alt)
            if ch is None:
                bad_ctx += 1
                continue
            counts[ch] += 1
        matrix[name] = counts
        total = sum(counts.values())
        meta[name] = {
            'snvs_considered': n,
            'dropped_unliftable': dropped_lift,
            'context_unavailable': bad_ctx,
            'assigned_to_channel': total,
            'ref_base_match_pct': round(100 * ref_ok / n, 4) if n else None,
        }
        print(f"  {name:<15} SNVs={n:>7,}  assigned={total:>7,}"
              f"  REF match={meta[name]['ref_base_match_pct']}%"
              + (f"  dropped(unliftable)={dropped_lift:,}" if lo else ''), flush=True)

    suffix = '_liftfiltered' if args.exclude_unliftable else ''
    with open(os.path.join(OUT, f'matrix96{suffix}.json'), 'w') as fh:
        json.dump({'matrix': matrix, 'meta': meta}, fh, indent=2)

    cols = [s for s in ORDER if s in matrix]
    with open(os.path.join(OUT, f'matrix96{suffix}.tsv'), 'w') as fh:
        fh.write('MutationType\t' + '\t'.join(cols) + '\n')
        for ch in CHANNELS:
            fh.write(ch + '\t' + '\t'.join(str(matrix[s][ch]) for s in cols) + '\n')

    print(f'\nWrote results/signatures/matrix96{suffix}.tsv  (96 x {len(cols)})')

    print('\n=== collapsed to 6 classes (% of assigned) ===')
    hdr = f"{'sample':<16}{'n':>9}" + ''.join(f'{s:>9}' for s in SUBS)
    print(hdr); print('-' * len(hdr))
    for s in cols:
        tot = sum(matrix[s].values()) or 1
        row = ''
        for sub in SUBS:
            v = sum(c for ch, c in matrix[s].items() if ch[2:5] == sub)
            row += f'{100 * v / tot:>9.2f}'
        print(f'{s:<16}{tot:>9,}' + row)

    print('\n=== most frequent channels per sample (top 4) ===')
    for s in cols:
        tot = sum(matrix[s].values()) or 1
        top = sorted(matrix[s].items(), key=lambda kv: -kv[1])[:4]
        print(f"  {s:<16}" + '  '.join(f'{ch} {100*v/tot:.1f}%' for ch, v in top))


if __name__ == '__main__':
    main()
