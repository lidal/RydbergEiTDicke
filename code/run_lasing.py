"""Analysis + figure for the U(1) (lasing) route and the three-transition comparison."""
import sys, os, json, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from lasing_route import *

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'figure.dpi': 200,
    'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5,
    'lines.linewidth': 1.4, 'savefig.bbox': 'tight', 'axes.axisbelow': True})
C = ['#1b4965', '#c1666b', '#4d908e', '#e09f3e', '#5f0f40']

hbar, e, a0 = 1.054571817e-34, 1.602176634e-19, 5.29177210903e-11
K = 1774.82 * e * a0 / hbar
MHz = 2 * np.pi * 1e6
kappa, Gamma = 2.0 * MHz, 0.02 * MHz      # bad cavity
out = {}

fig, ax = plt.subplots(1, 4, figsize=(13.6, 2.6))
fig.subplots_adjust(wspace=0.42)

# ---- (a) sensing curve: output power vs microwave field --------------------
w_cav, delta = 0.0, -3.0 * MHz            # bias so that E_mu tunes D through zero
E = np.linspace(0, 4.0, 700)              # mV/cm
Om = K * (E * 1e-3 * 1e2)
D = w_cav - (delta + Om / 2)
for i, r in enumerate([1.5, 3.0, 8.0]):
    G2 = r * kappa * Gamma
    I = np.array([steady_state(G2, kappa, Gamma, dd, 1, 1, 1)[0] for dd in D])
    ax[0].plot(E, I, color=plt.cm.viridis(i / 3.4), label=rf'$G^2/\kappa\Gamma={r:g}$')
ax[0].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[0].set_ylabel(r'output  $I/I_{\rm sat}$')
ax[0].legend(frameon=False)
ax[0].set_title('(a) lasing sensing curve', fontsize=8.5, loc='left')

# ---- (b) cavity-noise rejection -------------------------------------------
ratio = np.logspace(-1.3, 2.7, 300)       # kappa/Gamma
c_at = ratio / (ratio + 1)
c_cav = 1 / (ratio + 1)
ax[1].loglog(ratio, c_cav, color=C[1], label=r'$\partial\nu/\partial\omega_{\rm cav}$  (leak)')
ax[1].loglog(ratio, c_at, color=C[0], label=r'$\partial\nu/\partial\tilde\omega_0$  (signal)')
ax[1].axvline(kappa / Gamma, color='0.4', ls='--', lw=0.9)
ax[1].text(kappa / Gamma * 1.25, 2e-3, 'operating\npoint', fontsize=6.5, color='0.35')
ax[1].set_xlabel(r'$\kappa/\Gamma$  (bad-cavity parameter)')
ax[1].set_ylabel('pulling coefficient')
ax[1].legend(frameon=False, loc='center left')
ax[1].set_title('(b) cavity-drift immunity', fontsize=8.5, loc='left')
out['rejection_at_operating_point'] = float((kappa + Gamma) / Gamma)

# ---- (c) Rydberg interactions -> bistability (the NEPT, from the lasing side)
G2 = 3.0 * kappa * Gamma
D0g = np.linspace(0.0, 9.0, 260) * MHz
for i, chi in enumerate([0.0, 1.0, 2.0, 4.0]):
    lo, hi, xs = [], [], []
    for d0 in D0g:
        r = rydberg_selfconsistent(G2, kappa, Gamma, d0, chi * MHz)
        if r:
            xs.append(d0 / MHz); lo.append(min(r)); hi.append(max(r))
    if xs:
        ax[2].plot(xs, hi, color=plt.cm.plasma(i / 4.4), lw=1.4,
                   label=rf'$\chi/2\pi={chi:g}$ MHz')
        if chi > 0:
            ax[2].plot(xs, lo, color=plt.cm.plasma(i / 4.4), lw=0.9, ls=':')
ax[2].set_xlabel(r'bare detuning $D_0/2\pi$  (MHz)')
ax[2].set_ylabel(r'$I/I_{\rm sat}$')
ax[2].legend(frameon=False, fontsize=6.5)
ax[2].set_title('(c) interactions $\\rightarrow$ bistability', fontsize=8.5, loc='left')

# ---- (d) three-route sensitivity floors -----------------------------------
Ns = np.logspace(3, 8, 60)
g_common = Gamma                                    # matched decoherence
u1 = 2 * np.sqrt(2 * g_common / Ns) / K * 1e-2 * 1e9          # U(1) lasing
z2 = 2 * (0.72 * np.sqrt(g_common / (0.01 * Ns))) / K * 1e-2 * 1e9   # Z2, n_b = 0.01 N
z2_ideal = 2 * (0.72 * np.sqrt(g_common / Ns)) / K * 1e-2 * 1e9      # Z2 if n_b = N
ax[3].loglog(Ns, u1, color=C[0], label=r'$U(1)$ lasing  $\sqrt{2\Gamma/N}$')
ax[3].loglog(Ns, z2, color=C[1], label=r'$\mathbb{Z}_2$ Dicke  ($n_b=10^{-2}N$)')
ax[3].loglog(Ns, z2_ideal, color=C[1], ls=':', lw=1.0, label=r'$\mathbb{Z}_2$ if $n_b=N$')
ax[3].axhline(49, color=C[3], ls='--', lw=1.0)
ax[3].text(2e3, 60, 'Ding 2022 (NEPT)', fontsize=6.5, color=C[3])
ax[3].axhline(2.6, color=C[4], ls='--', lw=1.0)
ax[3].text(2e3, 3.2, 'Wang 2025 (NEPT + cavity)', fontsize=6.5, color=C[4])
ax[3].set_xlabel(r'atoms $N$')
ax[3].set_ylabel(r'$\sqrt{S_E}$  (nV cm$^{-1}$Hz$^{-1/2}$)')
ax[3].legend(frameon=False, fontsize=6.5)
ax[3].set_title(r'(d) floors at matched $\Gamma$', fontsize=8.5, loc='left')

fig.savefig('figs/fig5_lasing.pdf'); plt.close(fig)
print('fig5 done')

# ---------------------------------------------------------------- numbers
print("\n=== U(1) lasing route ===")
ca, cc = pulling_coefficients(kappa, Gamma)
print(f"  kappa/2pi={kappa/MHz} MHz, Gamma/2pi={Gamma/MHz} MHz  ->  kappa/Gamma = {kappa/Gamma:.0f}")
print(f"  pulling: signal {ca:.4f}, cavity leak {cc:.4e}  -> cavity drift rejected {1/cc:.0f}x")
print(f"  scale factor  d(nu)/d(E_mu) = {ca*(K/2)/(2*np.pi*1e6)*0.1:.4f} MHz per (mV/cm)"
      f"   [= {ca:.4f} x the half-AT slope {(K/2)/(2*np.pi*1e6)*0.1:.4f} MHz/(mV/cm)]")
rows = []
for N in [1e4, 1e5, 1e6, 1e7, 1e8]:
    dw, n, sE = sensitivity_pulling(kappa, Gamma, N, Gamma, K)
    z2v = 2*(0.72*np.sqrt(Gamma/(0.01*N)))/K*1e-2*1e9
    rows.append(dict(N=N, n_ph=n, linewidth_mHz=dw/2/np.pi*1e3, sE_nV=sE*1e9, z2_nV=z2v))
    print(f"   N={N:.0e}: n={n:.3e}, linewidth={dw/2/np.pi*1e3:9.3f} mHz, "
          f"S_E={sE*1e9:8.4f} nV/cm/rtHz   (Z2 at n_b=0.01N: {z2v:7.4f})")
out['rows'] = rows
out['pulling'] = dict(signal=ca, leak=cc)
json.dump(out, open('data/lasing.json','w'), indent=1)
print("\nwrote data/lasing.json")
