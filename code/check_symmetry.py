"""What breaks the |+> / |-> symmetry?  Nothing in the Hamiltonian -- the
reference state does.  Both branches couple to the cavity equally (lam/sqrt2).
The asymmetry enters only through which level the atoms actually sit in."""
import sys; sys.path.insert(0,'/home/user/RydbergEiTDicke/code')
import numpy as np
from two_branch import *

w, lam = 1.0*MHz, 0.35*MHz

def V_from_g(w, energies, g):
    """Potential matrix: cavity + one bosonic mode per cavity-coupled excitation."""
    n = len(energies)
    V = np.zeros((n+1, n+1))
    V[0,0] = w
    for i,e in enumerate(energies):
        V[i+1,i+1] = e
        V[0,i+1] = V[i+1,0] = 2*g
    return V

print("=== 1. The symmetric point delta = 0 ===")
print("   w_pm = +-Omega_mu/2 : the branches ARE symmetric about |g>.")
print("   Is |g> a minimum there?\n")
print("   Omega_mu/2pi   w_+/2pi   w_-/2pi   min eig V   |g> a minimum?")
for E in [0.5, 1.0, 2.0, 5.0]:
    O = Omega_mu(E); wp, wm = O/2, -O/2
    mn = np.min(np.linalg.eigvalsh(potential_matrix(w, wp, wm, lam)))
    print(f"    {O/MHz:9.4f}   {wp/MHz:+7.4f}   {wm/MHz:+7.4f}   {mn/MHz:+9.4f}       "
          f"{'yes' if mn>0 else 'NO'}")
print("\n   => at the symmetric point |g> sits BETWEEN the branches: a saddle, never")
print("      a minimum.  So exactly where + and - are symmetric there is NO transition.")

print("\n=== 2. Reference |g> (needs delta > Omega_mu/2: both branches above |g>) ===")
E = 2.0; O = Omega_mu(E); u = O/2
dp = ( 4*lam**2 + np.sqrt(16*lam**4 + 4*w**2*u**2) )/(2*w)
print(f"   E = {E} mV/cm: threshold at delta_+/2pi = {dp/MHz:.4f} MHz,"
      f"  w_-/2pi = {(dp-u)/MHz:.4f} MHz")
print(f"   asymptotically w_- -> 2 lam^2/w = {2*lam**2/w/MHz:.4f} MHz")

print("\n=== 3. The MIRROR case: reference |-> (delta < -Omega_mu/2, both below |g>) ===")
print("   Then only ONE excitation is cavity-coupled: |-> -> |g>, energy |w_-|,")
print("   coupling lam/sqrt2.  (|-> -> |+> is the MW transition, not cavity-coupled.)")
print("   That is a clean single-mode Dicke model, threshold  2 lam^2 = w |w_-|.\n")
print("   |w_-|/2pi   det V'      threshold |w_-|/2pi = 2 lam^2/w")
g = lam/np.sqrt(2)
for x in [0.10, 0.2450, 0.40, 1.0]:
    wm_abs = x*MHz
    Vp = np.array([[w, 2*g],[2*g, wm_abs]])
    print(f"   {x:9.4f}   {np.linalg.det(Vp)/MHz**2:+9.5f}   {2*lam**2/w/MHz:.4f}")
print("\n   => a genuine mirror transition exists, at the SAME |w_-| = 2 lam^2/w.")
print("      Its reference state is 'all atoms in |->' -- a Rydberg state.")
print("      That requires inversion, which is exactly the Route I obstacle.")
