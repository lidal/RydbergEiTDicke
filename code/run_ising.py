"""Dicke-Ising extension (Zhang et al. PRL 110, 090402): effect of Rydberg repulsion."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from dicke_core import mf_dicke_ising_2sub

w = 1.0
lams = np.linspace(0.05, 1.50, 146)

print("Cut at w0 = +1 (sensor operating point), z = 4:\n")
cuts = {}
for V in [0.0, 0.5, 2.0]:
    a = np.array([mf_dicke_ising_2sub(w, 1.0, L, V)['alpha'] for L in lams])
    onset = lams[int(np.argmax(a > 1e-4))]
    # first order <=> a finite jump that a sqrt onset cannot produce
    dj = np.max(np.diff(a))
    sqrt_ref = np.sqrt(lams[1] - lams[0]) * 1.0
    cuts[str(V)] = dict(onset=float(onset), max_jump=float(dj), alpha=a.tolist())
    print(f"  V = {V:4.1f}:  onset lam = {onset:.4f}  (lam_c^(V=0) = 0.5000),  "
          f"max step in |alpha| = {dj:.4f}")

print("\nPhase diagram in (w0, lam) at V = 2.0, z = 4:")
w0s = np.linspace(-2.0, 2.0, 41)
lam2 = np.linspace(0.05, 1.50, 41)
grid, labels = [], []
for ww in w0s:
    row = [mf_dicke_ising_2sub(w, ww, L, 2.0) for L in lam2]
    grid.append([r['alpha'] for r in row])
    labels.append([r['phase'] for r in row])
found = sorted({p for r in labels for p in r})
print(f"  phases present: {found}")
for ww, row in zip(w0s[::8], labels[::8]):
    seq = []
    for p in row:
        if not seq or seq[-1] != p:
            seq.append(p)
    print(f"   w0 = {ww:+5.2f}:  {' -> '.join(seq)}")

# first-order check via hysteresis: follow the minimum adiabatically up and down
print("\nHysteresis scan at w0 = -0.5, V = 2.0 (up vs down sweep in lam):")
def sweep(lamlist, w0v, V):
    from scipy.optimize import minimize
    out, t = [], np.array([0.1, 3.0])
    for L in lamlist:
        def en(tt):
            tA, tB = tt
            al = -(L / (2 * w)) * (np.sin(tA) + np.sin(tB))
            nA, nB = 0.5 * (1 - np.cos(tA)), 0.5 * (1 - np.cos(tB))
            return (w * al**2 - 0.25 * w0v * (np.cos(tA) + np.cos(tB))
                    + L * al * (np.sin(tA) + np.sin(tB)) + 0.5 * V * 4 * nA * nB)
        r = minimize(en, t, method='L-BFGS-B', bounds=[(0, np.pi)] * 2, tol=1e-14)
        t = r.x
        out.append(abs(-(L / (2 * w)) * (np.sin(t[0]) + np.sin(t[1]))))
    return np.array(out)

ls = np.linspace(0.05, 1.5, 200)
up, dn = sweep(ls, -0.5, 2.0), sweep(ls[::-1], -0.5, 2.0)[::-1]
hyst = float(np.max(np.abs(up - dn)))
print(f"   max |alpha_up - alpha_down| = {hyst:.5f}  -> "
      f"{'HYSTERETIC (first order)' if hyst > 1e-3 else 'no hysteresis (second order)'}")

json.dump(dict(lams=lams.tolist(), cuts=cuts, w0s=w0s.tolist(),
               lam2=lam2.tolist(), grid=grid, labels=labels,
               hysteresis=dict(lam=ls.tolist(), up=up.tolist(), down=dn.tolist(),
                               max_gap=hyst)),
          open('data/ising.json', 'w'), indent=1)
print("\nwrote data/ising.json")
