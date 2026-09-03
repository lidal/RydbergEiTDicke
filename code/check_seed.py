"""
Why the counter-rotating terms are the mechanism, not a refinement.

With the co-rotating (Tavis-Cummings) coupling alone the excitation number is
conserved, so |g,0> -- every atom in the ground state, no cavity photons -- is an
exact eigenstate and the cavity stays dark forever.  Restoring the counter-rotating
terms makes |g,0> non-stationary, and above threshold unstable.
"""
import sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'code')
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from dicke_core import spin_ops, boson_ops

def build(N, nmax, w, w0, lam, cr):
    """cr=0 -> Tavis-Cummings (co-rotating only); cr=1 -> full Dicke."""
    a, num, db = boson_ops(nmax)
    Jz, Jx, Jp, Jm, ds = spin_ops(N)
    Ib, Is = sp.identity(db,format='csr'), sp.identity(ds,format='csr')
    co = sp.kron(a, Jp) + sp.kron(a.T, Jm)          # a J+ + a^dag J-
    ct = sp.kron(a.T, Jp) + sp.kron(a, Jm)          # a^dag J+ + a J-   (counter-rotating)
    H = w*sp.kron(num,Is) + w0*sp.kron(Ib,Jz) + (lam/np.sqrt(N))*(co + cr*ct)
    return H.tocsr(), db, ds

N, nmax, w, w0 = 8, 30, 1.0, 1.0
lam = 1.2*0.5*np.sqrt(w*w0)          # 20% above the Dicke threshold

# |g,0> : all spins down (last spin basis state), zero photons
def psi0(db, ds):
    v = np.zeros(db*ds); v[0*ds + (ds-1)] = 1.0; return v

print("Is |g,0> an eigenstate?   residual || (H - <H>) |g,0> ||")
for cr, name in [(0,"Tavis-Cummings (coupling laser only)"), (1,"full Dicke (+ counter-rotating)")]:
    H, db, ds = build(N, nmax, w, w0, lam, cr)
    v = psi0(db, ds)
    Hv = H @ v
    E = v @ Hv
    res = np.linalg.norm(Hv - E*v)
    print(f"   {name:42s}  E = {E:+.4f},  residual = {res:.3e}")

print("\nTime evolution of <a^dag a> starting from |g,0>:")
from scipy.sparse.linalg import expm_multiply
for cr, name in [(0,"Tavis-Cummings"), (1,"full Dicke")]:
    H, db, ds = build(N, nmax, w, w0, lam, cr)
    a, num, _ = boson_ops(nmax); Is = sp.identity(ds,format='csr')
    Nop = sp.kron(num, Is).tocsr()
    v = psi0(db, ds).astype(complex)
    out=[]
    for t in [0.0, 1.0, 3.0, 10.0, 30.0]:
        vt = expm_multiply(-1j*H*t, v)
        out.append(np.real(np.vdot(vt, Nop@vt)))
    print(f"   {name:16s} t=0,1,3,10,30 -> " + "  ".join(f"{x:.3e}" for x in out))
