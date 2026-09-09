#!/usr/bin/env python3
"""
verify_call_sets.py -- data-integrity audit of Azenta somatic call sets.

Reports, per sample: total records, PASS records, PASS composition
(SNV / insertion / deletion / MNV), Ts/Tv on PASS and on all records,
and the count of records carrying the `germline` FILTER tag.

Also prints the tumour/normal pairing declared in each VCF header, which
is the authoritative record of the analysis design.

Usage:  python3 scripts/qc/verify_call_sets.py [data_dir]
"""
import gzip, glob, os, sys, collections

TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}


def header_design(path):
    """Extract tumour/normal sample names and pipeline provenance."""
    info = {'tumour': None, 'normal': None, 'sentieon': None,
            'reference': None, 'pon': False, 'samples': []}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            if line.startswith('##tumor_sample='):
                info['tumour'] = line.strip().split('=', 1)[1]
            elif line.startswith('##normal_sample='):
                info['normal'] = line.strip().split('=', 1)[1]
            elif line.startswith('##reference='):
                info['reference'] = line.strip().split('=', 1)[1]
            elif 'TNhaplotyper2' in line:
                if 'Version="' in line:
                    info['sentieon'] = line.split('Version="')[1].split('"')[0]
                if '--pon' in line:
                    info['pon'] = True
            elif line.startswith('#CHROM'):
                info['samples'] = line.rstrip('\n').split('\t')[9:]
    return info


def audit(path):
    counts = dict(total=0, passing=0, snv=0, ins=0, dele=0, mnv=0,
                  germline_tag=0, ts=0, tv=0, ts_all=0, tv_all=0)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.split('\t', 8)
            ref, alt, filt = cols[3], cols[4].split(',')[0], cols[6]
            counts['total'] += 1
            if 'germline' in filt:
                counts['germline_tag'] += 1
            is_snv = len(ref) == 1 and len(alt) == 1
            if is_snv:
                key = 'ts_all' if (ref, alt) in TRANSITIONS else 'tv_all'
                counts[key] += 1
            if filt != 'PASS':
                continue
            counts['passing'] += 1
            if is_snv:
                counts['snv'] += 1
                key = 'ts' if (ref, alt) in TRANSITIONS else 'tv'
                counts[key] += 1
            elif len(ref) == len(alt):
                counts['mnv'] += 1
            elif len(alt) > len(ref):
                counts['ins'] += 1
            else:
                counts['dele'] += 1
    return counts


def main(data_dir='data'):
    pattern = os.path.join(data_dir, '*', '*', '*_somatic.vcf.gz')
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f'No somatic VCFs found under {data_dir}')

    print('\n=== ANALYSIS DESIGN (from VCF headers) ===')
    for path in files:
        d = header_design(path)
        pairing = (f"paired: tumour={d['tumour']} normal={d['normal']}"
                   if d['normal'] else f"tumour-only: {d['tumour']}")
        print(f"  {os.path.basename(path)}")
        print(f"    {pairing}")
        print(f"    sentieon={d['sentieon']}  PoN={'yes' if d['pon'] else 'NO'}"
              f"  samples={d['samples']}")

    print('\n=== CALL SET COMPOSITION ===')
    hdr = (f"{'sample':<16}{'total':>11}{'PASS':>9}{'%PASS':>7}{'SNV':>9}"
           f"{'INS':>8}{'DEL':>8}{'MNV':>7}{'Ts/Tv':>8}{'Ts/Tv(all)':>12}"
           f"{'germline':>11}")
    print(hdr)
    print('-' * len(hdr))
    for path in files:
        name = os.path.basename(path).replace('_somatic.vcf.gz', '')
        c = audit(path)
        ratio = lambda a, b: f'{a / b:.2f}' if b else 'NA'
        print(f"{name:<16}{c['total']:>11,}{c['passing']:>9,}"
              f"{100 * c['passing'] / c['total']:>6.1f}%{c['snv']:>9,}"
              f"{c['ins']:>8,}{c['dele']:>8,}{c['mnv']:>7,}"
              f"{ratio(c['ts'], c['tv']):>8}{ratio(c['ts_all'], c['tv_all']):>12}"
              f"{c['germline_tag']:>11,}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'data')
