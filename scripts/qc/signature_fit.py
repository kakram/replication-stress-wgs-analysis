#!/usr/bin/env python3
"""
signature_fit.py -- fit the 96-channel matrix against COSMIC SBS signatures.

Method
------
Non-negative least squares against the COSMIC v3.4 GRCh38 reference
signatures, followed by sparsity pruning: signatures contributing less than
MIN_CONTRIB of the total are dropped and the fit repeated, which is the
standard remedy for the over-fitting that occurs when ~80 signatures are
fitted to a few thousand mutations. Reconstruction quality is reported as
cosine similarity between observed and fitted spectra; below about 0.90 a
fit should not be trusted.

Because these call sets are small (roughly 6,000-74,000 substitutions),
contributions are bootstrapped (resampling mutations with replacement) and
reported with 95% intervals. Wide intervals are expected and should be
quoted rather than hidden.

Also reports, independently of any fit, the cosine similarity between each
sample's raw spectrum and each individual COSMIC signature. That is a
model-free summary and is more robust than a decomposition at this sample
size.

Caveat: COSMIC signatures were derived from primary tumours. These are
long-passaged immortalised cell lines and culture-specific processes are
not represented in the catalogue, so an imperfect fit is expected.
"""
import json, os, sys

import numpy as np
from scipy.optimize import nnls

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(REPO, 'results', 'signatures')
COSMIC = ('/sessions/rcw-01utqydzxgsezjbb66mkbisx/.local/lib/python3.10/'
          'site-packages/SigProfilerAssignment/data/Reference_Signatures/'
          'GRCh38/COSMIC_v3.4_SBS_GRCh38.txt')

MIN_CONTRIB = 0.05
N_BOOT = 200
RNG = np.random.default_rng(20260904)
ORDER = ['WTUN', 'WTAPH', 'B4UN', 'B4APH',
         'WT-U_cleaned', 'WT-A_cleaned', 'G3-U_cleaned', 'G3-A_cleaned']


def load_cosmic():
    with open(COSMIC) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        names = header[1:]
        rows, channels = [], []
        for line in fh:
            p = line.rstrip('\n').split('\t')
            channels.append(p[0])
            rows.append([float(x) for x in p[1:]])
    return np.array(rows), names, channels


def cosine(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na and nb else 0.0


def sparse_fit(obs, S, min_contrib=MIN_CONTRIB):
    """NNLS then prune weak signatures and refit. Returns (weights, active)."""
    active = np.arange(S.shape[1])
    for _ in range(20):
        w, _ = nnls(S[:, active], obs)
        total = w.sum()
        if total <= 0:
            return np.zeros(len(active)), active
        keep = w / total >= min_contrib
        if keep.all():
            return w, active
        if not keep.any():
            best = int(np.argmax(w))
            active = active[[best]]
            continue
        active = active[keep]
    w, _ = nnls(S[:, active], obs)
    return w, active


def main():
    S, names, channels = load_cosmic()
    with open(os.path.join(OUT, 'matrix96.json')) as fh:
        blob = json.load(fh)
    matrix = blob['matrix']

    out = {'cosmic_version': 'COSMIC v3.4 (GRCh38)', 'samples': {}}
    print(f'COSMIC v3.4 GRCh38: {S.shape[1]} signatures x {S.shape[0]} channels\n')

    for name in [s for s in ORDER if s in matrix]:
        obs = np.array([matrix[name][c] for c in channels], float)
        n = obs.sum()
        w, active = sparse_fit(obs, S)
        recon = S[:, active] @ w
        cs = cosine(obs, recon)
        contrib = {names[i]: float(w[j] / w.sum())
                   for j, i in enumerate(active) if w.sum() > 0}

        boots = {k: [] for k in contrib}
        probs = obs / n
        for _ in range(N_BOOT):
            draw = RNG.multinomial(int(n), probs).astype(float)
            bw, ba = sparse_fit(draw, S)
            tot = bw.sum()
            bmap = {names[i]: bw[j] / tot for j, i in enumerate(ba)} if tot else {}
            for k in boots:
                boots[k].append(bmap.get(k, 0.0))
        ci = {k: (round(100 * float(np.percentile(v, 2.5)), 1),
                  round(100 * float(np.percentile(v, 97.5)), 1))
              for k, v in boots.items()}

        sims = sorted(((cosine(obs, S[:, i]), names[i]) for i in range(S.shape[1])),
                      reverse=True)[:5]

        out['samples'][name] = {
            'n_snv': int(n),
            'cosine_similarity_of_fit': round(cs, 4),
            'contributions_pct': {k: round(100 * v, 2)
                                  for k, v in sorted(contrib.items(),
                                                     key=lambda kv: -kv[1])},
            'contributions_ci95': ci,
            'closest_single_signatures': [(s, round(c, 4)) for c, s in sims],
        }

        print(f'{name}  (n={int(n):,}, reconstruction cosine={cs:.3f})')
        for k, v in sorted(contrib.items(), key=lambda kv: -kv[1]):
            lo, hi = ci[k]
            print(f'    {k:<8} {100*v:>5.1f}%   95% CI {lo:>5.1f} - {hi:<5.1f}')
        print('    closest single signatures: ' +
              ', '.join(f'{s} {c:.3f}' for c, s in sims[:3]))
        print()

    with open(os.path.join(OUT, 'signature_fit.json'), 'w') as fh:
        json.dump(out, fh, indent=2)
    print('Wrote results/signatures/signature_fit.json')


if __name__ == '__main__':
    main()
