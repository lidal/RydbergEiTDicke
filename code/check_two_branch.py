"""
Exact two-branch threshold and the true Protocol B scale factor.

Raman couples |g> to the BARE |r> = (|+> + |->)/sqrt2, hence to both dressed
branches with lambda/sqrt2 each.  Normal-phase stability:

    2 lam^2 (1/w_+ + 1/w_-) = w      with  w_pm = delta +- u,  u = Omega_mu/2
 => 4 lam^2 delta = w (delta^2 - u^2) = w w_+ w_-                      (threshold)

Differentiating at fixed lam, w:   ddelta/du = w u / (w delta - 2 lam^2).
The single-branch model instead gives ddelta/du = -1 (on |+>) or +1 (on |->).
"""
import numpy as np

MHz = 2*np.pi*1e6
K = 1774.82*1.602176634e-19*5.29177210903e-11/1.054571817e-34
def Om(E): return K*(E*1e-3*1e2)
def E_of(O): return O/K*1e-2*1e3

w, lam = 1.0*MHz, 0.35*MHz

def delta_th(u):
    """Solve w delta^2 - 4 lam^2 delta - w u^2 = 0, take the root with w_- > 0."""
    a, b, c = w, -4*lam**2, -w*u**2
    disc = b*b - 4*a*c
    r = [(-b+np.sqrt(disc))/(2*a), (-b-np.sqrt(disc))/(2*a)]
    good = [x for x in r if x - u > 0]
    return max(good) if good else np.nan

def slope(u):
    d = delta_th(u)
    return w*u/(w*d - 2*lam**2)

print("Exact two-branch threshold, w/2pi = 1 MHz, lam/2pi = 0.35 MHz")
print("(2 lam^2/w /2pi = %.4f MHz is the |-> -dominated threshold splitting)\n"
      % (2*lam**2/w/MHz))
print("  E_mu[mV/cm]  Om_mu/2pi   delta/2pi   w_+/2pi   w_-/2pi   ddelta_th/dOm_mu   vs +1/2")
for E in [0.05, 0.1, 0.22, 0.5, 1.0, 2.0, 5.0, 20.0]:
    u = Om(E)/2
    d = delta_th(u)
    s = slope(u)/2.0                       # d delta_th / d Omega_mu
    print(f"   {E:7.2f}    {2*u/MHz:8.4f}   {d/MHz:8.4f}  {(d+u)/MHz:8.4f}  {(d-u)/MHz:8.4f}"
          f"      {s:+8.5f}        {s/0.5:6.3f}")

print("\nScale factor referred to field, ddelta_th/dE_mu  [MHz per (mV/cm)]:")
print("   single-branch model in the paper:  -1.1355  (built on |+>)")
for E in [0.22, 1.0, 5.0, 20.0]:
    u = Om(E)/2
    s = slope(u)/2.0*K*1e-2*1e3            # rad/s per (mV/cm)
    print(f"   two-branch, E = {E:5.2f} mV/cm :  {s/MHz:+8.4f}")

print("\nlam_c comparison at fixed delta (does the threshold itself move?):")
for E in [0.22, 1.0, 5.0]:
    u = Om(E)/2; d = delta_th(u)
    l2 = np.sqrt(w*(d+u)*(d-u)/(4*d))      # from 4 lam^2 delta = w w_+ w_-
    l1p = 0.5*np.sqrt(w*(d+u))             # single-branch on |+>
    l1m = 0.5*np.sqrt(w*(d-u))             # single-branch on |->
    print(f"   E={E:5.2f}: lam_c(2br)={l2/MHz:.4f}  lam_c(1br,|+>)={l1p/MHz:.4f}  "
          f"lam_c(1br,|->)={l1m/MHz:.4f} MHz")
