"""Fig 1, corrected: both MW-dressed branches retained."""
import sys, os, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from two_branch import *

plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':7.5,
 'xtick.labelsize':8,'ytick.labelsize':8,'figure.dpi':200,'axes.grid':True,
 'grid.alpha':0.25,'grid.linewidth':0.5,'lines.linewidth':1.5,'savefig.bbox':'tight',
 'axes.axisbelow':True})
C=['#1b4965','#c1666b','#4d908e','#e09f3e']
w, lam, delta0 = 1.0*MHz, 0.35*MHz, 3.0*MHz

fig, ax = plt.subplots(1,3, figsize=(10.4,2.7)); fig.subplots_adjust(wspace=0.40)

# (a) phase boundary in (E_mu, lambda) at fixed delta
E = np.linspace(0, 5.5, 700)
lc = np.array([lam_c(w, delta0, Omega_mu(x))/MHz for x in E])
ax[0].plot(E, lc, color=C[0])
ax[0].fill_between(E, lc, 1.0, color=C[1], alpha=0.10)
ax[0].fill_between(E, 0, lc, color=C[0], alpha=0.08)
ax[0].axhline(lam/MHz, color='0.35', ls='--', lw=1.0)
ax[0].text(3.4,0.62,'superradiant\n(bright)',fontsize=7,color=C[1],ha='center')
ax[0].text(1.1,0.14,'normal (dark)',fontsize=7,color=C[0],ha='center')
ax[0].set_xlabel(r'microwave field $E_\mu$  (mV/cm)')
ax[0].set_ylabel(r'$\lambda/2\pi$  (MHz)'); ax[0].set_ylim(0,0.9); ax[0].set_xlim(0,5.5)
ax[0].set_title(r'(a) boundary at fixed $\delta$: field turns it ON', fontsize=8.5, loc='left')

# (b) Protocol B lock point vs field, with the linear asymptote
Eg = np.linspace(0.02, 6.0, 500)
dth = np.array([delta_th(w, lam, Omega_mu(x))/MHz for x in Eg])
ax[1].plot(Eg, dth, color=C[0], label='two branches (correct)')
asym = 2*lam**2/w/MHz + 0.5*K/(2*np.pi*1e6)*0.1*Eg
ax[1].plot(Eg, asym, color=C[3], ls='--', lw=1.1, label=r'asymptote, slope $+d/2\hbar$')
ax[1].plot(Eg, 0.49 - 1.1355*Eg, color=C[1], ls=':', lw=1.2,
           label='single branch (earlier, wrong)')
ax[1].set_xlabel(r'$E_\mu$  (mV/cm)'); ax[1].set_ylabel(r'$\delta_{\rm th}/2\pi$  (MHz)')
ax[1].legend(frameon=False, fontsize=6.6); ax[1].set_ylim(-1.5, 7)
ax[1].set_title('(b) lock point: sign was wrong', fontsize=8.5, loc='left')

# (c) scale factor approaches d/2hbar only for large field
sc = np.array([slope_dth_dOm(w, lam, Omega_mu(x))*K/(2*np.pi*1e6)*0.1 for x in Eg])
ax[2].plot(Eg, sc, color=C[0])
ax[2].axhline(1.1355, color=C[3], ls='--', lw=1.1)
ax[2].text(3.0, 1.02, r'$d/2\hbar = 1.1355$', fontsize=7, color=C[3])
for xv, lab in [(0.5,'92%'), (1.0,'98%'), (2.0,'99.4%')]:
    yv = slope_dth_dOm(w,lam,Omega_mu(xv))*K/(2*np.pi*1e6)*0.1
    ax[2].plot([xv],[yv],'o',color=C[1],ms=4)
    ax[2].annotate(lab,(xv,yv),(xv+0.25,yv-0.22),fontsize=6.5,color=C[1])
ax[2].set_xlabel(r'$E_\mu$  (mV/cm)')
ax[2].set_ylabel(r'$\partial\delta_{\rm th}/\partial E_\mu$  [MHz/(mV/cm)]')
ax[2].set_ylim(0,1.3)
ax[2].set_title(r'(c) linear only for $\Omega_\mu\gg2\lambda^2/\omega$', fontsize=8.5, loc='left')

fig.savefig('figs/fig1_threshold.pdf'); fig.savefig('figs/fig1_threshold.png', dpi=185)
plt.close(fig); print('fig1 regenerated (two-branch)')
