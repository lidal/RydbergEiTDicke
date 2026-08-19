"""Gain, bandwidth and input-referred sensitivity of the critical-Dicke electrometer."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
from open_dicke import (drift, lam_crit, steady_state, response_and_noise, slow_rate)

# ---------------------------------------------------------------- constants
hbar = 1.054571817e-34
e    = 1.602176634e-19
a0   = 5.29177210903e-11
d_ea0 = 1774.82                      # 87Rb 53D5/2 -> 54P3/2
d    = d_ea0 * e * a0                # C m
d_over_hbar = d / hbar               # (rad/s) per (V/m)
print(f"dipole d = {d_ea0} e a0 = {d:.4e} C m")
print(f"Omega_mu/E = {d_over_hbar:.4e} rad/s per V/m "
      f"= {d_over_hbar/2/np.pi/1e6:.4f} MHz/(V/m) "
      f"= {d_over_hbar/2/np.pi/1e6*0.1:.4f} MHz per (mV/cm)")

# --------------------------------------------------- experimental scale set
MHz = 2 * np.pi * 1e6
w      = 1.0 * MHz          # effective cavity detuning in the Raman frame
kappa  = 0.05 * MHz         # cavity HWHM  (kappa/2pi = 50 kHz)
gamma  = 0.005 * MHz        # collective-mode decoherence (gamma/2pi = 5 kHz)
w0     = 1.0 * MHz          # bias splitting  w0 = delta + Omega_LO/2
drive  = 1.0 * MHz          # cavity drive strength (sets <b>)

lc = lam_crit(w, w0, kappa, gamma)
print(f"\nlam_c/2pi = {lc/2/np.pi/1e6:.6f} MHz  (undamped {0.5*np.sqrt(w*w0)/2/np.pi/1e6:.6f} MHz)")

# ------------------------------------- gain, bandwidth, input-referred noise
eps = np.logspace(-4, -0.3, 60)        # eps = 1 - lam/lam_c
rows = []
for ee in eps:
    lam = lc * (1 - ee)
    vbar = steady_state(w, w0, lam, kappa, gamma, drive)
    nb = 0.5 * (vbar[2] ** 2 + vbar[3] ** 2)      # <b^dag b> of the collective mode

    # optimise the homodyne angle for input-referred noise
    best = None
    for phi in np.linspace(0, np.pi, 181):
        T, So, Sw = response_and_noise(w, w0, lam, kappa, gamma, drive, phi, nu=0.0)
        if best is None or Sw < best[2]:
            best = (T, So, Sw, phi)
    T, So, Sw, phi = best
    bw = slow_rate(w, w0, lam, kappa, gamma)

    # gain: dc response of the measured output quadrature per unit dw0,
    # normalised to unit drive
    rows.append(dict(eps=float(ee), lam=float(lam), T=float(T), S_out=float(So),
                     S_w0=float(Sw), bandwidth=float(bw), phi=float(phi),
                     nb=float(nb), gain_bw=float(T * bw)))

print("\n   eps        |T| (gain)     bandwidth/2pi[kHz]   G*BW        S_out/vac   sqrt(S_w0) [Hz/rtHz]")
for r in rows[::8]:
    print(f"  {r['eps']:.2e}   {r['T']:.4e}   {r['bandwidth']/2/np.pi/1e3:10.4f}   "
          f"{r['gain_bw']:.4e}   {r['S_out']/0.5:8.4f}   {np.sqrt(r['S_w0'])/2/np.pi:.4e}")

gb = np.array([r['gain_bw'] for r in rows])
print(f"\ngain x bandwidth: mean {gb.mean():.4e}, spread {(gb.max()-gb.min())/gb.mean()*100:.2f}%")

# ---------------------------------- how the input-referred noise scales with nb
print("\ninput-referred noise vs collective-mode occupation (eps = 1e-3):")
lam = lc * (1 - 1e-3)
for dr in [0.1, 0.3, 1.0, 3.0, 10.0]:
    vb = steady_state(w, w0, lam, kappa, gamma, dr * MHz)
    nb = 0.5 * (vb[2] ** 2 + vb[3] ** 2)
    best = min(response_and_noise(w, w0, lam, kappa, gamma, dr * MHz, p, 0.0)[2]
               for p in np.linspace(0, np.pi, 181))
    print(f"   drive={dr:5.1f} MHz  <b^dag b>={nb:11.3f}  "
          f"sqrt(S_w0)/2pi = {np.sqrt(best)/2/np.pi:.4e} Hz/rtHz   "
          f"product*sqrt(nb) = {np.sqrt(best*nb)/2/np.pi:.4e}")

# ------------------------------------------------- squeezing of the output
print("\noutput-quadrature noise (min over phi) vs eps  [vacuum = 1]:")
for ee in [1e-1, 1e-2, 1e-3, 1e-4]:
    lam = lc * (1 - ee)
    sm = min(response_and_noise(w, w0, lam, kappa, gamma, drive, p, 0.0)[1]
             for p in np.linspace(0, np.pi, 361))
    print(f"   eps={ee:.0e}   S_out/S_vac = {sm/0.5:.4f}  ({10*np.log10(sm/0.5):+.2f} dB)")

json.dump(dict(w=w, w0=w0, kappa=kappa, gamma=gamma, drive=drive, lam_c=lc,
               d_ea0=d_ea0, d_over_hbar=d_over_hbar, rows=rows),
          open('data/sensing.json', 'w'), indent=1)
print("\nwrote data/sensing.json")
