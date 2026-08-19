"""
The decisive test: hold the collective-mode occupation n_b fixed (this is the
real physical resource, bounded by n_b << N) and ask whether approaching the
critical point still improves the input-referred sensitivity.
"""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.optimize import brentq
from open_dicke import lam_crit, steady_state, response_and_noise, slow_rate

MHz = 2 * np.pi * 1e6
w = w0 = 1.0 * MHz
kappa, gamma = 0.05 * MHz, 0.005 * MHz
lc = lam_crit(w, w0, kappa, gamma)

def nb_of(lam, drive):
    v = steady_state(w, w0, lam, kappa, gamma, drive)
    return 0.5 * (v[2] ** 2 + v[3] ** 2)

def drive_for_nb(lam, target):
    f = lambda ld: np.log(nb_of(lam, np.exp(ld))) - np.log(target)
    return np.exp(brentq(f, np.log(1e-12 * MHz), np.log(1e9 * MHz), xtol=1e-12))

def best_Sw0(lam, drive, nu=0.0, n_th=0.0):
    return min(response_and_noise(w, w0, lam, kappa, gamma, drive, p, nu, n_th)[2]
               for p in np.linspace(0, np.pi, 361))

TARGET = 1.0e4          # fixed collective-mode occupation
print(f"Fixed n_b = {TARGET:.0e};  lam_c/2pi = {lc/MHz:.6f} MHz\n")
print("   eps       lam/lam_c    drive/2pi[MHz]     sqrt(S_w0)/2pi [Hz/rtHz]   bandwidth/2pi[kHz]")
res = []
for ee in [3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5]:
    lam = lc * (1 - ee)
    dr = drive_for_nb(lam, TARGET)
    S = best_Sw0(lam, dr)
    bw = slow_rate(w, w0, lam, kappa, gamma)
    res.append(dict(eps=ee, S=float(np.sqrt(S) / 2 / np.pi), bw=float(bw / MHz * 1e3)))
    print(f"  {ee:.1e}    {1-ee:.6f}   {dr/MHz:14.6e}   {np.sqrt(S)/2/np.pi:20.6f}   {bw/MHz*1e3:12.4f}")

s = np.array([r['S'] for r in res])
print(f"\n  -> input-referred noise varies by only {(s.max()-s.min())/s.mean()*100:.3f}% "
      f"over 4.5 decades in eps.")

# same test at finite Fourier frequency (a real receiver has an IF)
print("\nAt finite IF (nu/2pi = 10 kHz), fixed n_b:")
print("   eps        sqrt(S_w0)/2pi [Hz/rtHz]")
for ee in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
    lam = lc * (1 - ee)
    dr = drive_for_nb(lam, TARGET)
    S = best_Sw0(lam, dr, nu=2 * np.pi * 10e3)
    print(f"  {ee:.1e}    {np.sqrt(S)/2/np.pi:14.6f}")

# scaling of the floor with n_b and with the atomic decoherence rate
print("\nFloor vs n_b (eps = 1e-3):")
lam = lc * (1 - 1e-3)
for nb in [1e2, 1e3, 1e4, 1e5, 1e6]:
    dr = drive_for_nb(lam, nb)
    S = np.sqrt(best_Sw0(lam, dr)) / 2 / np.pi
    print(f"   n_b={nb:8.0e}   sqrt(S_w0)/2pi = {S:12.6f} Hz/rtHz    x sqrt(n_b) = {S*np.sqrt(nb):.6f}")

print("\nFloor vs atomic decoherence gamma (eps = 1e-3, n_b = 1e4):")
for gg in [0.0005, 0.005, 0.05, 0.5]:
    gam = gg * MHz
    lcg = lam_crit(w, w0, kappa, gam)
    lam = lcg * (1 - 1e-3)
    f = lambda ld: np.log(0.5 * np.sum(steady_state(w, w0, lam, kappa, gam, np.exp(ld))[2:] ** 2)) - np.log(1e4)
    dr = np.exp(brentq(f, np.log(1e-12 * MHz), np.log(1e9 * MHz), xtol=1e-12))
    S = min(response_and_noise(w, w0, lam, kappa, gam, dr, p, 0.0)[2]
            for p in np.linspace(0, np.pi, 361))
    print(f"   gamma/2pi={gg*1e3:7.2f} kHz   sqrt(S_w0)/2pi = {np.sqrt(S)/2/np.pi:12.6f} Hz/rtHz"
          f"    / sqrt(gamma/2pi) = {np.sqrt(S)/2/np.pi/np.sqrt(gg*1e6):.6f}")
