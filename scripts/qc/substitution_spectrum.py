#!/usr/bin/env python3
"""
substitution_spectrum.py -- pyrimidine-normalised substitution spectrum and
an orientation-bias (8-oxoG) diagnostic for the filtered call sets.

Why
---
Ts/Tv in these call sets sits below 1.0 and does not rise under stricter
filtering. Two explanations compete: (a) a genuine transversion-heavy
mutational process, or (b) a residual technical artefact -- classically
8-oxoguanine damage during library preparation, which produces C>A (G>T)
transversions.

The two are separable. Oxidative artefacts are strand-asymmetric: the
alternate allele is observed overwhelmingly on one read orientation. Real
mutations are present on both. TNhaplotyper2 records F1R2 and F2R1 read
counts per allele, so the asymmetry can be measured directly.

An orientation ratio near 0.50 means balanced support (consistent with real
variants). A ratio approaching 0 or 1 indicates orientation bias.

Filters applied match tier T4 of filter_cascade.py, and T5 where paired.
"""
import glob, gzip, json, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(REPO, 'data')
OUT = os.path.join(REPO, 'results', 'cascade')

MIN_DP, MIN_VAF, MIN_ALT, MIN_POPAF = 20, 0.05, 3, 3.0
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
CLASSES = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
POPAF_RE = re.compile(r'(?:^|;)POPAF=([^;]+)')
ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']


def pyrimidine(ref, alt):
    if ref in 'CT':
        return f'{ref}>{alt}'
    return f'{COMP[ref]}>{COMP[alt]}'


def run(vcf_path):
    name = os.path.basename(vcf_path).replace('_somatic.vcf.gz', '')
    spec = {c: 0 for c in CLASSES}
    orient = {c: [0, 0] for c in CLASSES}      # [F1R2 alt, F2R1 alt]
    paired = None

    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                paired = len(line.rstrip('\n').split('\t')) > 10
                continue
            c = line.rstrip('\n').split('\t')
            ref, alt = c[3], c[4].split(',')[0]
            if c[6] != 'PASS' or len(ref) != 1 or len(alt) != 1:
                continue
            if ref not in COMP or alt not in COMP:
                continue
            keys = c[8].split(':')
            t = dict(zip(keys, c[9].split(':')))
            try:
                if int(t.get('DP', 0)) < MIN_DP:
                    continue
                ad = t['AD'].split(',')
                if int(ad[1]) < MIN_ALT:
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
                    if int(n['AD'].split(',')[1]) != 0:   # T5: strictly private
                        continue
                except (KeyError, ValueError, IndexError):
                    continue
            m = POPAF_RE.search(c[7])
            if m and float(m.group(1).split(',')[0]) < MIN_POPAF:
                continue

            k = pyrimidine(ref, alt)
            spec[k] += 1
            try:
                orient[k][0] += int(t['F1R2'].split(',')[1])
                orient[k][1] += int(t['F2R1'].split(',')[1])
            except (KeyError, ValueError, IndexError):
                pass

    total = sum(spec.values()) or 1
    rec = {
        'sample': name, 'paired': paired, 'n_snv': sum(spec.values()),
        'spectrum': spec,
        'spectrum_pct': {k: round(100 * v / total, 2) for k, v in spec.items()},
        'orientation_f1r2_fraction': {
            k: (round(v[0] / (v[0] + v[1]), 4) if (v[0] + v[1]) else None)
            for k, v in orient.items()},
    }
    with open(os.path.join(OUT, f'{name}_spectrum.json'), 'w') as fh:
        json.dump(rec, fh, indent=2)
    return rec


if __name__ == '__main__':
    recs = [run(p) for p in
            sorted(glob.glob(os.path.join(DATA, '*', '*', '*_somatic.vcf.gz')))]
    recs.sort(key=lambda r: ORDER.index(r['sample'])
              if r['sample'] in ORDER else 99)

    print('=== SUBSTITUTION SPECTRUM (% of filtered SNVs) ===')
    h = f"{'sample':<15}{'n':>9}" + ''.join(f'{c:>9}' for c in CLASSES)
    print(h); print('-' * len(h))
    for r in recs:
        print(f"{r['sample']:<15}{r['n_snv']:>9,}" +
              ''.join(f"{r['spectrum_pct'][c]:>9}" for c in CLASSES))

    print('\n=== ORIENTATION BIAS: fraction of alt reads in F1R2 ===')
    print('    (0.50 = balanced / real;  ->0 or ->1 = strand artefact)')
    print(h.replace(f"{'n':>9}", ' ' * 9)); print('-' * len(h))
    for r in recs:
        row = ''
        for c in CLASSES:
            v = r['orientation_f1r2_fraction'][c]
            row += f"{v if v is not None else '-':>9}"
        print(f"{r['sample']:<15}{'':>9}" + row)
