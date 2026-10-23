import os
import shutil
import sys

import numpy as np

if len(sys.argv) < 3:
    sys.exit(f"usage: {sys.argv[0]} <keys> <out_dir> [h=0.45]")
keys_file, root = sys.argv[1], sys.argv[2]
h = float(sys.argv[3]) if len(sys.argv) > 3 else 0.45

K, MU_D, NQ, L = 32, 0.45, 1_000_000, 1 << 10
SEED_D, SEED_Q = 202, 303

keys = np.fromfile(keys_file, dtype=np.uint64)[1:]
N = len(keys)
cuts = (np.arange(K + 1) * N) // K
D = np.random.default_rng(SEED_D).permutation(np.linspace(MU_D - h, MU_D + h, K))
rng = np.random.default_rng(SEED_Q)
counts = np.floor(NQ * np.diff(cuts) / N).astype(np.int64)
counts[:NQ - counts.sum()] += 1


def is_empty(x):
    j = np.searchsorted(keys, x, side='left')
    return (j == N) | (keys[np.minimum(j, N - 1)] > x + np.uint64(L - 1))


qs = []
for g in range(K):
    lo, hi, need, got = cuts[g], cuts[g + 1], counts[g], []
    tries = 0
    while sum(len(a) for a in got) < need and tries < 200:
        m = max(4 * (need - sum(len(a) for a in got)), 1000)
        k = keys[rng.integers(lo, hi, size=m)]
        span = np.uint64(min(1 << 62, int(round(2 ** (30 * (1 - D[g]))))))
        off = rng.integers(0, span, size=m, dtype=np.uint64, endpoint=True)
        x = k + off
        ok = (x >= k) & is_empty(x)
        got.append(x[ok])
        tries += 1
    x = np.concatenate(got)[:need]
    qs.append(x)
    print(f"region {g:2d}: keys {hi - lo}, D {D[g]:.2f}, queries {len(x)}")

q = np.concatenate(qs)
np.random.default_rng(SEED_Q + 1).shuffle(q)


def with_header(a, dtype):
    return np.concatenate([np.array([len(a)], dtype=dtype), a.astype(dtype)])


for r in (0, 5, 10):
    out = os.path.join(root, f"r{r}")
    os.makedirs(out, exist_ok=True)
    shutil.copy(keys_file, os.path.join(out, "keys"))
    with_header(q, np.uint64).tofile(os.path.join(out, "left"))
    with_header(q + np.uint64((1 << r) - 1), np.uint64).tofile(os.path.join(out, "right"))
    with_header(np.zeros(len(q)), np.int32).tofile(os.path.join(out, "result"))
print(f"wrote {len(q)} queries to {root}")
