"""SciPy's MINRES on the exported fixtures, at a fixed step count with the
residual test switched off, so what is compared is the recurrence.

    python dev_notes/spikes/minres_scipy.py

rtol = 0 removes the only tolerance MINRES exposes. SciPy still stops on its own
breakdown and conditioning tests, which is the point on the two fixtures named
for those, so the iteration count it actually ran is written out beside the
iterate.
"""
import os
import numpy as np
from scipy.sparse.linalg import minres

d = os.environ.get("MINRES_XDIR",
                   os.path.join(os.environ.get("TMPDIR", "/tmp"), "minres-xcheck"))
print("dir:", d)

cases = open(os.path.join(d, "cases.txt")).read().split()
steps = [int(s) for s in open(os.path.join(d, "steps.txt")).read().split()]

for name in cases:
    A = np.loadtxt(os.path.join(d, name + "-A.txt"))
    b = np.loadtxt(os.path.join(d, name + "-b.txt"))
    for k in steps:
        count = [0]

        def bump(_xk, c=count):
            c[0] += 1

        x, info = minres(A, b, rtol=0.0, maxiter=k, callback=bump)
        np.savetxt(os.path.join(d, "%s-scipy-%d.txt" % (name, k)), x)
        with open(os.path.join(d, "%s-scipyit-%d.txt" % (name, k)), "w") as f:
            f.write("%d\n" % count[0])
    print(name, "done", A.shape)
