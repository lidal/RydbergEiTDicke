"""Finite-size scaling of gap, photon number and QFI at the Dicke critical point."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from dicke_core import gap, photon_number, qfi_w0, ground_state

w, w0 = 1.0, 1.0
lam_c = 0.5 * np.sqrt(w * w0)

# --- nmax convergence check at the largest N we will use -------------------
print("nmax convergence (N=200, lam=lam_c):")
for nm in [40, 60, 80, 120]:
    q = qfi_w0(200, nm, w, w0, lam_c)
    n, _ = photon_number(200, nm, w, w0, lam_c)
    print(f"   nmax={nm:4d}  n={n:.6f}  QFI={q:.6f}")

Ns = [4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]
rows = []
for N in Ns:
    nm = 100
    g, _ = gap(N, nm, w, w0, lam_c)
    n, jz = photon_number(N, nm, w, w0, lam_c)
    q = qfi_w0(N, nm, w, w0, lam_c)
    rows.append(dict(N=N, gap=float(g), n=float(n), jz_per_N=float(jz / N), qfi=float(q)))
    print(f"N={N:4d}  gap={g:.6f}  n={n:.6f}  QFI={q:.6f}")

A = np.array([[r['N'], r['gap'], r['n'], r['qfi']] for r in rows])
sel = A[:, 0] >= 24          # fit the asymptotic tail
fit = {}
for name, col in [('gap', 1), ('n', 2), ('qfi', 3)]:
    p = np.polyfit(np.log(A[sel, 0]), np.log(A[sel, col]), 1)
    fit[name] = dict(exponent=float(p[0]), prefactor=float(np.exp(p[1])))
    print(f"fit  {name:4s}:  exponent = {p[0]:+.4f}  (expected "
          f"{ {'gap': -1/3, 'n': 1/3, 'qfi': 4/3}[name]:+.4f})")

# --- QFI vs coupling for several N (peak position / height) ----------------
lams = np.linspace(0.20, 0.80, 121)
curves = {}
for N in [16, 32, 64, 128, 256]:
    qs = [qfi_w0(N, 100, w, w0, L) for L in lams]
    curves[str(N)] = qs
    i = int(np.argmax(qs))
    print(f"N={N:4d}  QFI peak {qs[i]:.3f} at lam/lam_c={lams[i]/lam_c:.4f}")

# --- QFI vs distance-from-threshold in the normal phase (thermodynamic-limit
#     prediction F_Q = (1/2) eps^-2 (d eps/d w0)^2, eps = 1 - lam^2/lam_c^2 ) --
eps_grid = np.logspace(-2.2, -0.2, 25)
normal = []
for e in eps_grid:
    L = lam_c * np.sqrt(1 - e)
    q = qfi_w0(512, 140, w, w0, L)
    deps_dw0 = 4 * L ** 2 / (w * w0 ** 2)
    pred = 0.5 * deps_dw0 ** 2 / e ** 2
    normal.append(dict(eps=float(e), qfi=float(q), gaussian_pred=float(pred)))

out = dict(w=w, w0=w0, lam_c=lam_c, scaling=rows, fits=fit,
           lams=lams.tolist(), qfi_vs_lam=curves, normal_phase=normal)
os.makedirs('data', exist_ok=True)
json.dump(out, open('data/scaling.json', 'w'), indent=1)
print("wrote data/scaling.json")
