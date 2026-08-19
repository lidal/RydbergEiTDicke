"""Realistic sensitivity estimates and the Dicke-Ising (Rydberg-interaction) extension."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from scipy.optimize import brentq
from open_dicke import lam_crit, steady_state, response_and_noise, slow_rate
from dicke_core import mf_dicke_ising, mf_order_parameter

hbar, e, a0 = 1.054571817e-34, 1.602176634e-19, 5.29177210903e-11
d_ea0 = 1774.82                                   # 87Rb 53D5/2 -> 54P3/2
d = d_ea0 * e * a0
K = d / hbar                                      # rad/s per V/m
MHz = 2 * np.pi * 1e6

w = w0 = 1.0 * MHz
kappa, gamma = 0.05 * MHz, 0.005 * MHz
lc = lam_crit(w, w0, kappa, gamma)

def nb_of(lam, dr):
    v = steady_state(w, w0, lam, kappa, gamma, dr)
    return 0.5 * (v[2] ** 2 + v[3] ** 2)

def drive_for_nb(lam, tgt):
    return np.exp(brentq(lambda ld: np.log(nb_of(lam, np.exp(ld))) - np.log(tgt),
                         np.log(1e-12 * MHz), np.log(1e12 * MHz), xtol=1e-13))

def sens(lam, dr, nu=0.0):
    S = min(response_and_noise(w, w0, lam, kappa, gamma, dr, p, nu)[2]
            for p in np.linspace(0, np.pi, 361))
    sw = np.sqrt(S)                      # rad/s per rtHz on w0
    sE = 2 * sw / K                      # V/m per rtHz   (dw0 = dOmega/2, E = hbar Omega/d)
    return sw, sE * 1e-2                 # V/cm per rtHz

print("Sensitivity vs distance from threshold at fixed collective occupation n_b = 1e4")
print("  eps      lam/lam_c   sqrt(S_w0)/2pi [Hz/rtHz]   S_E [nV/cm/rtHz]   BW/2pi [kHz]")
best = None
table = []
for ee in [0.9, 0.5, 0.3, 0.1, 3e-2, 1e-2, 5e-3, 3e-3, 2e-3, 1e-3, 3e-4, 1e-4]:
    lam = lc * (1 - ee)
    dr = drive_for_nb(lam, 1e4)
    sw, sE = sens(lam, dr)
    bw = slow_rate(w, w0, lam, kappa, gamma) / MHz * 1e3
    table.append(dict(eps=ee, S_w0_Hz=sw / 2 / np.pi, S_E_nV=sE * 1e9, bw_kHz=bw))
    print(f"  {ee:.1e}   {1-ee:.6f}   {sw/2/np.pi:16.4f}   {sE*1e9:16.4f}   {bw:10.3f}")
    if best is None or sE < best[1]:
        best = (ee, sE, bw)
print(f"\n  optimum: eps = {best[0]:.1e}  ->  {best[1]*1e9:.4f} nV/cm/rtHz, "
      f"bandwidth {best[2]:.2f} kHz")
ref = table[0]['S_E_nV']
print(f"  improvement vs far-from-threshold (eps=0.9): {ref/ (best[1]*1e9):.1f}x")

print("\nSensitivity vs the number of atoms in the collective mode (eps = 3e-3):")
lam = lc * (1 - 3e-3)
nb_rows = []
for nb in [1e3, 1e4, 1e5, 1e6, 1e7]:
    dr = drive_for_nb(lam, nb)
    sw, sE = sens(lam, dr)
    sql = np.sqrt(gamma / nb)                       # N indep. atoms, T2 = 1/gamma
    nb_rows.append(dict(nb=nb, S_E_nV=sE*1e9, ratio_to_sql=float(sw/sql)))
    print(f"   n_b = {nb:8.0e}   S_E = {sE*1e9:10.4f} nV/cm/rtHz   "
          f"(atomic SQL for N=n_b would be {2*sql/K*1e-2*1e9:8.4f});  ratio {sw/sql:.3f}")

# ------------------------------------------------ null / threshold read-out
print("\nThreshold (null) electrometry:  lam_c = sqrt(w (delta + Omega_mu/2))/2")
print("   Omega_mu/2pi [MHz]   E [mV/cm]   lam_c/2pi [MHz]   d lam_c/dE [kHz per (uV/cm)]")
thr = []
for E_mVcm in [0.0, 0.01, 0.1, 1.0]:
    Om = K * (E_mVcm * 1e-3 * 1e2)                  # mV/cm -> V/m
    ww0 = 0.0 + Om / 2
    lcth = 0.5 * np.sqrt(w * ww0) if ww0 > 0 else 0.0
    dl = 0.25 * np.sqrt(w / ww0) * K / 2 * 1e-4 if ww0 > 0 else np.inf   # per uV/cm
    thr.append(dict(E_mVcm=E_mVcm, Om_MHz=Om / MHz, lam_c_MHz=lcth / MHz))
    print(f"   {Om/MHz:16.5f}   {E_mVcm:9.3f}   {lcth/MHz:14.5f}   "
          f"{dl/2/np.pi/1e3 if np.isfinite(dl) else float('inf'):12.4f}")

# --------------------------------- Rydberg interactions: Dicke-Ising, 1st order
print("\nDicke-Ising mean field (infinite-range Rydberg repulsion V):")
wI, w0I = 1.0, 1.0
lam_scan = np.linspace(0.30, 0.90, 4001)
di = {}
for V in [0.0, 0.5, 2.0, 5.0]:
    nb_arr = np.array([mf_dicke_ising(wI, w0I, L, V)[2] for L in lam_scan])
    jump = np.max(np.diff(nb_arr))
    ith = int(np.argmax(nb_arr > 1e-9))
    di[str(V)] = dict(lam_th=float(lam_scan[ith]), max_jump=float(jump),
                      first_order=bool(jump > 1e-3))
    print(f"   V = {V:4.1f}:  threshold lam/sqrt(w w0) = {lam_scan[ith]:.4f}, "
          f"largest jump in n_bar = {jump:.5f}  -> "
          f"{'FIRST order' if jump > 1e-3 else 'second order'}")

json.dump(dict(d_ea0=d_ea0, K=K, lam_c=lc, w=w, w0=w0, kappa=kappa, gamma=gamma,
               eps_table=table, nb_table=nb_rows, threshold=thr, dicke_ising=di,
               best=dict(eps=best[0], S_E_nV=best[1]*1e9, bw_kHz=best[2])),
          open('data/numbers.json', 'w'), indent=1)
print("\nwrote data/numbers.json")
