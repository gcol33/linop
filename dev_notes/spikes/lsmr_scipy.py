"""SciPy's LSMR on the exported fixtures, at a fixed step count with every
stopping test switched off, so what is compared is the recurrence.

    python dev_notes/spikes/lsmr_scipy.py
"""
import os
import numpy as np
from scipy.sparse.linalg import lsmr

d = os.environ.get("LSMR_XDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "lsmr-xcheck"))
print("dir:", d)

cases = open(os.path.join(d, "cases.txt")).read().split()
steps = [int(s) for s in open(os.path.join(d, "steps.txt")).read().split()]

for name in cases:
    A = np.loadtxt(os.path.join(d, name + "-A.txt"))
    b = np.loadtxt(os.path.join(d, name + "-b.txt"))
    for k in steps:
        out = lsmr(A, b, damp=0.0, atol=0.0, btol=0.0, conlim=0.0, maxiter=k)
        x = out[0]
        np.savetxt(os.path.join(d, "%s-scipy-%d.txt" % (name, k)), x)
    print(name, "done", A.shape)
