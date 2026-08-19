"""Figures for the critical-Dicke Rydberg electrometer note."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from dicke_core import mf_order_parameter
from open_dicke import lam_crit, steady_state, response_and_noise, slow_rate

plt.rcParams.update({
    'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'figure.dpi': 200,
    'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5,
    'lines.linewidth': 1.4, 'savefig.bbox': 'tight', 'axes.axisbelow': True,
})
C = ['#1b4965', '#c1666b', '#4d908e', '#e09f3e', '#5f0f40']
os.makedirs('figs', exist_ok=True)

hbar, e, a0 = 1.054571817e-34, 1.602176634e-19, 5.29177210903e-11
K = 1774.82 * e * a0 / hbar
MHz = 2 * np.pi * 1e6

# =====================================================================  FIG 1
# Sensing curve: the critical coupling is a function of the microwave field.
fig, ax = plt.subplots(1, 2, figsize=(7.2, 2.6))
fig.subplots_adjust(wspace=0.42)
w = 1.0
E = np.linspace(0, 1.5, 400)                     # mV/cm
Om = K * (E * 1e-3 * 1e2)                        # rad/s
w0_phys = Om / 2
lc_phys = 0.5 * np.sqrt(w * MHz * w0_phys)
ax[0].plot(E, lc_phys / MHz, color=C[0])
ax[0].fill_between(E, 0, lc_phys / MHz, color=C[0], alpha=0.10)
ax[0].fill_between(E, lc_phys / MHz, 0.9, color=C[1], alpha=0.10)
ax[0].text(1.0, 0.13, 'normal\n(dark cavity)', ha='center', fontsize=7.5, color=C[0])
ax[0].text(0.35, 0.62, 'superradiant\n(bright cavity)', ha='center', fontsize=7.5, color=C[1])
ax[0].axhline(0.35, color='0.35', ls='--', lw=1.0)
ax[0].plot([0.365], [0.35], 'o', color=C[3], ms=5, zorder=5)
ax[0].annotate('operating point', (0.365, 0.35), (0.60, 0.30),
               fontsize=7, arrowprops=dict(arrowstyle='->', lw=0.7, color='0.35'))
ax[0].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[0].set_ylabel(r'collective coupling $\lambda/2\pi$  (MHz)')
ax[0].set_xlim(0, 1.5); ax[0].set_ylim(0, 0.9)
ax[0].set_title(r'(a)  threshold line $\lambda_c=\frac{1}{2}\sqrt{\omega\,\Omega_\mu/2}$',
                fontsize=8.5, loc='left')

lam_fix = 0.35 * MHz
w0_grid = np.linspace(1e-4, 3.0, 800) * MHz
nb = np.array([mf_order_parameter(w * MHz, x, lam_fix)[1] for x in w0_grid])
E_grid = (2 * w0_grid) / K * 1e-2 * 1e3
ax[1].plot(E_grid, nb, color=C[1])
ax[1].fill_between(E_grid, 0, nb, color=C[1], alpha=0.15)
w0c = 4 * lam_fix ** 2 / (w * MHz)
Ec = (2 * w0c) / K * 1e-2 * 1e3
ax[1].axvline(Ec, color='0.35', ls='--', lw=1.0)
ax[1].text(Ec * 1.05, 0.10, r'$E_\mu^{c}$', fontsize=8, color='0.3')
ax[1].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[1].set_ylabel(r'$\bar n = \langle a^\dagger a\rangle/N$')
ax[1].set_xlim(0, 3 * Ec); ax[1].set_ylim(bottom=0)
ax[1].set_title(r'(b)  order parameter at fixed $\lambda$', fontsize=8.5, loc='left')
fig.savefig('figs/fig1_threshold.pdf'); plt.close(fig)
print('fig1 done')

# =====================================================================  FIG 2
# Finite-size scaling of gap and QFI at criticality.
S = json.load(open('data/scaling.json'))
N = np.array([r['N'] for r in S['scaling']])
gap = np.array([r['gap'] for r in S['scaling']])
qfi = np.array([r['qfi'] for r in S['scaling']])
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.35))
fig.subplots_adjust(wspace=0.38)
ax[0].loglog(N, gap, 'o', ms=3.5, color=C[0], label='exact diagonalisation')
ax[0].loglog(N, gap[-1] * (N / N[-1]) ** (-1 / 3), '-', color=C[1], lw=1.1,
             label=r'$\propto N^{-1/3}$')
ax[0].set_xlabel('$N$'); ax[0].set_ylabel(r'energy gap $\Delta/\omega$')
ax[0].legend(frameon=False); ax[0].set_title('(a) critical slowing down', fontsize=8.5, loc='left')

ax[1].loglog(N, qfi, 'o', ms=3.5, color=C[0], label='exact diagonalisation')
ax[1].loglog(N, qfi[-1] * (N / N[-1]) ** (4 / 3), '-', color=C[1], lw=1.1,
             label=r'$\propto N^{4/3}$')
ax[1].loglog(N, qfi[-1] * (N / N[-1]), '--', color='0.5', lw=1.0,
             label=r'$\propto N$  (SQL)')
ax[1].set_xlabel('$N$'); ax[1].set_ylabel(r'$F_Q[\tilde\omega_0]\ \ (\omega^{-2})$')
ax[1].legend(frameon=False); ax[1].set_title('(b) ground-state QFI at $\\lambda_c$', fontsize=8.5, loc='left')

lams = np.array(S['lams']); lc0 = S['lam_c']
for i, n in enumerate(['16', '32', '64', '128', '256']):
    ax[2].semilogy(lams / lc0, S['qfi_vs_lam'][n], color=plt.cm.viridis(i / 4.6),
                   lw=1.2, label=f'$N={n}$')
ax[2].axvline(1.0, color='0.4', ls='--', lw=0.9)
ax[2].set_xlabel(r'$\lambda/\lambda_c$'); ax[2].set_ylabel(r'$F_Q[\tilde\omega_0]$')
ax[2].legend(frameon=False, ncol=1); ax[2].set_xlim(0.45, 1.6)
ax[2].set_title('(c) QFI across the transition', fontsize=8.5, loc='left')
fig.savefig('figs/fig2_scaling.pdf'); plt.close(fig)
print('fig2 done')

# =====================================================================  FIG 3
# Open system: gain, bandwidth and input-referred sensitivity.
w = w0 = 1.0 * MHz
kappa, gamma = 0.05 * MHz, 0.005 * MHz
lc = lam_crit(w, w0, kappa, gamma)

def nb_of(lam, dr):
    v = steady_state(w, w0, lam, kappa, gamma, dr)
    return 0.5 * (v[2] ** 2 + v[3] ** 2)

def drive_for_nb(lam, tgt):
    return np.exp(brentq(lambda ld: np.log(nb_of(lam, np.exp(ld))) - np.log(tgt),
                         np.log(1e-14 * MHz), np.log(1e14 * MHz), xtol=1e-13))

eps = np.logspace(-5, -0.05, 70)
gain, bw, sE, sqz = [], [], [], []
for ee in eps:
    lam = lc * (1 - ee)
    dr = drive_for_nb(lam, 1e4)
    best = min((response_and_noise(w, w0, lam, kappa, gamma, dr, p, 0.0)[2], p)
               for p in np.linspace(0, np.pi, 361))
    sE.append(2 * np.sqrt(best[0]) / K * 1e-2 * 1e9)          # nV/cm/rtHz
    T = response_and_noise(w, w0, lam, kappa, gamma, MHz, np.pi / 2, 0.0)[0]
    gain.append(T)
    bw.append(slow_rate(w, w0, lam, kappa, gamma) / MHz * 1e3)
    smin = min(response_and_noise(w, w0, lam, kappa, gamma, dr, p, 0.0)[1]
               for p in np.linspace(0, np.pi, 361))
    sqz.append(10 * np.log10(smin / 0.5))

fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.35))
fig.subplots_adjust(wspace=0.38)
a0_ = ax[0]; a1 = a0_.twinx()
a0_.loglog(eps, gain, color=C[0]); a0_.set_ylabel(r'transduction gain $|T|$', color=C[0])
a0_.tick_params(axis='y', colors=C[0])
a1.loglog(eps, bw, color=C[1]); a1.set_ylabel('bandwidth (kHz)', color=C[1])
a1.tick_params(axis='y', colors=C[1]); a1.grid(False)
a0_.set_xlabel(r'$\varepsilon = 1-\lambda/\lambda_c$')
a0_.set_title('(a) gain vs bandwidth', fontsize=8.5, loc='left')

ax[1].loglog(eps, sE, color=C[0])
i = int(np.argmin(sE))
ax[1].plot(eps[i], sE[i], 'o', color=C[3], ms=5, zorder=5)
ax[1].annotate(f'{sE[i]:.2f} nV/cm/$\\sqrt{{\\rm Hz}}$\nat $\\varepsilon={eps[i]:.1e}$',
               (eps[i], sE[i]), (eps[i] * 3, sE[i] * 4), fontsize=6.5,
               arrowprops=dict(arrowstyle='->', lw=0.7, color='0.4'))
ax[1].set_xlabel(r'$\varepsilon = 1-\lambda/\lambda_c$')
ax[1].set_ylabel(r'$\sqrt{S_E}$  (nV cm$^{-1}$Hz$^{-1/2}$)')
ax[1].set_title(r'(c) sensitivity at fixed $n_b=10^4$', fontsize=8.5, loc='left')

ax[2].semilogx(eps, sqz, color=C[0])
ax[2].axhline(0, color='0.4', ls='--', lw=0.9)
ax[2].set_xlabel(r'$\varepsilon = 1-\lambda/\lambda_c$')
ax[2].set_ylabel('output noise (dB rel. vacuum)')
ax[2].set_title('(c) ponderomotive squeezing', fontsize=8.5, loc='left')
fig.savefig('figs/fig3_open.pdf'); plt.close(fig)
print('fig3 done  (best', min(sE), 'nV/cm/rtHz )')

json.dump(dict(eps=eps.tolist(), gain=gain, bw_kHz=bw, sE_nV=sE, sqz_dB=sqz),
          open('data/fig3.json', 'w'), indent=1)
