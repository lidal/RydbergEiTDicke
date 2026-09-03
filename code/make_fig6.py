"""Fig 6: what Protocol B looks like on the detector."""
import sys, os, json, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from protocol_b import *
from dicke_core import photon_number

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'figure.dpi': 200,
    'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5,
    'lines.linewidth': 1.5, 'savefig.bbox': 'tight', 'axes.axisbelow': True})
C = ['#1b4965', '#c1666b', '#4d908e', '#e09f3e']

w, lam = 1.0 * MHz, 0.35 * MHz
fig, ax = plt.subplots(1, 4, figsize=(13.8, 2.75))
fig.subplots_adjust(wspace=0.40)

# ---- (a) the read-out trace: cavity emission vs two-photon detuning --------
d = np.linspace(0.0, 0.70, 1400) * MHz
fields = [0.0, 0.15, 0.30]
for i, E in enumerate(fields):
    w0t = d + Omega_mu(E) / 2
    n = nbar_mf(w, w0t, lam)
    ax[0].plot(d / MHz, n, color=plt.cm.viridis(i / 3.3), label=f'$E_\\mu={E:.2f}$ mV/cm')
    dt = delta_threshold(w, lam, E) / MHz
    ax[0].plot([dt], [0], 'v', color=plt.cm.viridis(i / 3.3), ms=6, clip_on=False)
ax[0].annotate('lock here; the edge moves\n$-1.14$ MHz per (mV/cm)', (0.16, 0.012), (0.24, 0.075),
               fontsize=6.8, arrowprops=dict(arrowstyle='->', lw=0.7, color='0.4'))
ax[0].set_xlabel(r'two-photon detuning $\delta/2\pi$  (MHz)')
ax[0].set_ylabel(r'emitted photons per atom  $\bar n$')
ax[0].legend(frameon=False, loc='upper right')
ax[0].set_ylim(0, 0.135); ax[0].set_xlim(0, 0.70)
ax[0].set_title('(a) Protocol B read-out: cavity emission', fontsize=8.5, loc='left')

# ---- (b) how sharp is the edge? finite-N exact diagonalisation -------------
dz = np.linspace(0.36, 0.60, 46) * MHz
ax[1].plot(dz / MHz, nbar_mf(w, dz, lam), color='0.25', lw=1.8,
           label=r'mean field ($N\to\infty$)')
for i, N in enumerate([32, 128, 512]):
    y = [photon_number(N, 90, 1.0, max(x / MHz, 1e-6), 0.35)[0] / N for x in dz]
    ax[1].plot(dz / MHz, y, color=plt.cm.plasma(i / 3.6), lw=1.2, label=f'$N={N}$')
ax[1].axvline(4 * lam ** 2 / w / MHz, color=C[1], ls='--', lw=0.9)
ax[1].set_xlabel(r'two-photon detuning $\delta/2\pi$  (MHz)')
ax[1].set_ylabel(r'$\bar n$')
ax[1].legend(frameon=False, fontsize=6.8)
ax[1].set_title('(b) the edge is rounded by finite $N$', fontsize=8.5, loc='left')

# ---- (c) calibration line --------------------------------------------------
Es = np.linspace(0, 0.43, 100)
ax[2].plot(Es, delta_threshold(w, lam, Es) / MHz, color=C[0])
for E in fields:
    ax[2].plot([E], [delta_threshold(w, lam, E) / MHz], 'o', color=C[3], ms=5, zorder=5)
ax[2].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[2].set_ylabel(r'lock point $\delta_{\rm th}/2\pi$  (MHz)')
ax[2].text(0.05, 0.09, r'$\dfrac{\partial\delta_{\rm th}}{\partial E_\mu}=-\dfrac{d}{2\hbar}$'
                       '\nexact, dipole-referenced', fontsize=7.5, color=C[0])
ax[2].set_title('(c) calibration: a frequency, not a fit', fontsize=8.5, loc='left')

# ---- (d) what the conventional read-out would show at the same fields ------
g21, g31, g41 = 3.0 * MHz, 0.1 * MHz, 0.1 * MHz
Om_c = 5.0 * MHz
d2 = np.linspace(-1.6, 1.6, 1600) * MHz
OD = 3.0
bare = eit_at_transmission(np.array([50.0 * MHz]), Om_c, 0.0, g21, g31, g41)[0]
for i, E in enumerate(fields):
    a = eit_at_transmission(d2, Om_c, Omega_mu(E), g21, g31, g41)
    ax[3].plot(d2 / MHz, np.exp(-OD * a / bare), color=plt.cm.viridis(i / 3.3),
               label=f'$E_\\mu={E:.2f}$ mV/cm')
ax[3].set_xlabel(r'two-photon detuning $\delta/2\pi$  (MHz)')
ax[3].set_ylabel('probe transmission')
ax[3].legend(frameon=False, fontsize=6.8)
ax[3].set_title('(d) conventional EIT-AT, same fields', fontsize=8.5, loc='left')
ax[3].set_xlim(-1.6, 1.6)

os.makedirs('figs', exist_ok=True)
fig.savefig('figs/fig6_protocolB.pdf')
fig.savefig('figs/fig6_protocolB.png', dpi=190)
plt.close(fig)
print('fig6 done')

# ---- quantitative edge sharpness ------------------------------------------
# The mean-field edge is a KINK (nbar rises linearly below threshold), so it has
# no intrinsic width. Finite N rounds it over the critical region, whose width in
# the coupling is eps ~ N^{-2/3}; translated to detuning that is w0t_th * N^{-2/3}.
w0t_th = 4 * lam ** 2 / w
print("\nRounding of the edge (tail on the dark side of threshold):")
print("   N     nbar at delta_th    delta-delta_th where nbar halves    predicted w0t_th*N^(-2/3)")
for N in [32, 128, 512]:
    n_at = photon_number(N, 90, 1.0, w0t_th / MHz, 0.35)[0] / N
    xs = np.linspace(w0t_th, w0t_th + 0.25 * MHz, 120)
    ys = np.array([photon_number(N, 90, 1.0, x / MHz, 0.35)[0] / N for x in xs])
    k = np.where(ys < 0.5 * n_at)[0]
    meas = (xs[k[0]] - w0t_th) / MHz if len(k) else float('nan')
    print(f"  {N:5d}    {n_at:.5f}              {meas:8.4f} MHz                   "
          f"{w0t_th/MHz * N**(-2/3):8.4f} MHz")
print("\n   extrapolating the N^(-2/3) law: N=1e6 -> %.1f Hz, N=1e8 -> %.1f Hz"
      % (w0t_th/MHz*1e6**(-2/3)*1e6, w0t_th/MHz*1e8**(-2/3)*1e6))
print("   => quantum rounding is negligible for realistic N; the edge width is set")
print("      by the atomic decoherence gamma, not by finite-size effects.")

print("\nAT doublet splitting for the same fields (for comparison):")
for E in fields:
    print(f"   E={E:.2f} mV/cm -> AT splitting {Omega_mu(E)/MHz:.3f} MHz vs "
          f"EIT linewidth ~{(Om_c**2/(4*g21))/MHz:.2f} MHz"
          )
