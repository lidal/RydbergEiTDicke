"""Fig 6: what Protocol B looks like on the detector."""
import sys, os, json, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from protocol_b import eit_at_transmission
from two_branch import *
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
w, lam = 1.0*MHz, 0.35*MHz
d = np.linspace(0.35, 3.0, 2400)*MHz
fields = [1.0, 1.5, 2.0]
for i, E in enumerate(fields):
    O = Omega_mu(E)
    y = np.array([nbar_mf(w, x, O, lam) for x in d])
    ax[0].plot(d/MHz, y, color=plt.cm.viridis(i/3.3), label=f'$E_\\mu={E:.1f}$ mV/cm')
    dt = delta_th(w, lam, O)/MHz
    ax[0].plot([dt], [0], 'v', color=plt.cm.viridis(i/3.3), ms=6, clip_on=False)
ax[0].annotate('Protocol B locks\nto this edge', (2.53, 0.004), (1.55, 0.075),
               fontsize=6.6, arrowprops=dict(arrowstyle='->', lw=0.7, color='0.35'))
ax[0].annotate('', xy=(1.135, 0.128), xytext=(1.407, 0.128),
               arrowprops=dict(arrowstyle='<->', lw=0.8, color=C[1]))
ax[0].text(1.27, 0.133, r'$2\lambda^2/\omega$', fontsize=6.4, color=C[1], ha='center')
ax[0].set_xlabel(r'two-photon detuning $\delta/2\pi$  (MHz)')
ax[0].set_ylabel(r'emitted photons per atom  $\bar n$')
ax[0].legend(frameon=False, loc='center right', fontsize=6.8)
ax[0].set_xlim(0.35, 3.0); ax[0].set_ylim(0, 0.152)
ax[0].set_title('(a) a superradiant window that slides', fontsize=8.5, loc='left')

# ---- (b) how sharp is the edge? finite-N exact diagonalisation -------------
# Generic critical-scaling result: the edge is rounded over ~ w0_th * N^(-2/3).
# Computed in the single-branch Dicke model, where exact diagonalisation is cheap;
# the N^(-2/3) law is universal for this class and carries over unchanged.
def nbar_1br(w_, w0_, lam_):
    w0_ = np.asarray(w0_, float)
    mu = w_*w0_/(4*lam_**2)
    return np.where((mu < 1) & (w0_ > 0), (lam_**2/w_**2)*(1-np.clip(mu,0,1)**2), 0.0)

dz = np.linspace(0.36, 0.60, 46)*MHz
ax[1].plot(dz/MHz, nbar_1br(1.0*MHz, dz, 0.35*MHz), color='0.25', lw=1.8,
           label=r'mean field ($N\to\infty$)')
for i, N in enumerate([32, 128, 512]):
    y = [photon_number(N, 90, 1.0, max(x/MHz, 1e-6), 0.35)[0]/N for x in dz]
    ax[1].plot(dz/MHz, y, color=plt.cm.plasma(i/3.6), lw=1.2, label=f'$N={N}$')
ax[1].axvline(4*(0.35*MHz)**2/(1.0*MHz)/MHz, color=C[1], ls='--', lw=0.9)
ax[1].set_xlabel(r'detuning from threshold  (MHz)')
ax[1].set_ylabel(r'$\bar n$')
ax[1].legend(frameon=False, fontsize=6.8)
ax[1].set_title('(b) finite-$N$ rounding, $\\propto N^{-2/3}$', fontsize=8.5, loc='left')

# ---- (c) calibration line, corrected ---------------------------------------
Es = np.linspace(0.03, 3.0, 400)
dth = np.array([delta_th(w, lam, Omega_mu(x))/MHz for x in Es])
ax[2].plot(Es, dth, color=C[0], label='two branches (correct)')
ax[2].plot(Es, 2*lam**2/w/MHz + 1.1355*Es, color=C[3], ls='--', lw=1.1,
           label=r'asymptote $+d/2\hbar$')
ax[2].plot(Es, 0.49 - 1.1355*Es, color=C[1], ls=':', lw=1.2, label='single branch (wrong)')
ax[2].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[2].set_ylabel(r'lock point $\delta_{\rm th}/2\pi$  (MHz)')
ax[2].legend(frameon=False, fontsize=6.5)
ax[2].set_ylim(-1.2, 4.2)
ax[2].set_title('(c) calibration: sign corrected', fontsize=8.5, loc='left')

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

print("\nProtocol B scale factor vs field (asymptote 1.1355 MHz/(mV/cm)):")
for E in [0.1, 0.5, 1.0, 2.0]:
    sE = slope_dth_dOm(w, lam, Omega_mu(E))*K/(2*np.pi*1e6)*0.1
    print(f"   E={E:4.1f} mV/cm -> {sE:+.4f}  ({100*sE/1.1355:.1f}% of d/2hbar)")
