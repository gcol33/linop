"""SciPy's LSMR and LSQR run to convergence on the exported fixtures, so the
forward error each reaches can be compared against this implementation's.

    python dev_notes/spikes/lsmr_converged.py
"""
import os
import numpy as np
from scipy.sparse.linalg import lsmr, lsqr

d = os.environ.get("LSMR_XDIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "lsmr-xcheck"))
cases = open(os.path.join(d, "cases.txt")).read().split()

for name in cases:
    A = np.loadtxt(os.path.join(d, name + "-A.txt"))
    b = np.loadtxt(os.path.join(d, name + "-b.txt"))
    truth = np.loadtxt(os.path.join(d, name + "-truth.txt"))
    row = [name]
    for label in ("lsmr", "lsqr"):
        if label == "lsmr":
            out = lsmr(A, b, atol=1e-12, btol=1e-12, conlim=1e16, maxiter=5000)
        else:
            out = lsqr(A, b, atol=1e-12, btol=1e-12, conlim=1e16, iter_lim=5000)
        x = out[0]
        itn = out[2]
        err = np.max(np.abs(x - truth)) / np.max(np.abs(truth))
        r = b - A @ x
        bw = np.linalg.norm(A.T @ r) / (np.linalg.norm(A, 2) * np.linalg.norm(r))
        row.append("%s %4d it  fwd %.2e  bw %.2e" % (label, itn, err, bw))
        np.savetxt(os.path.join(d, "%s-scipy-conv-%s.txt" % (name, label)), x)
    print("  " + "  |  ".join(row))
