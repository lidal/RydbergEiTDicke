"""Fig 7: the two roots of the threshold equation, and why only one is a critical point."""
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

w, lam, E = 1.0*MHz, 0.35*MHz, 2.0
O = Omega_mu(E); u = O/2
a,b,c = w, -4*lam**2, -w*u**2
disc=b*b-4*a*c
dp=(-b+np.sqrt(disc))/(2*a); dm=(-b-np.sqrt(disc))/(2*a)

d = np.linspace(-5, 5, 3000)*MHz
det = np.array([np.linalg.det(potential_matrix(w, x+u, x-u, lam)) for x in d])
mn  = np.array([np.min(np.linalg.eigvalsh(potential_matrix(w, x+u, x-u, lam))) for x in d])

fig, ax = plt.subplots(1, 2, figsize=(8.6, 2.8)); fig.subplots_adjust(wspace=0.32)

ax[0].plot(d/MHz, det/MHz**3, color=C[0])
ax[0].axhline(0, color='0.5', lw=0.8)
for r, lab in [(dm,r'$\delta_-$'), (dp,r'$\delta_+$')]:
    ax[0].axvline(r/MHz, color=C[1], ls='--', lw=1.0)
    ax[0].text(r/MHz, 13, lab, color=C[1], fontsize=8, ha='center')
ax[0].set_xlabel(r'$\delta/2\pi$  (MHz)'); ax[0].set_ylabel(r'$\det V$  (MHz$^3$)')
ax[0].set_title(r'(a) two roots: $\delta_+\delta_-=-(\Omega_\mu/2)^2$', fontsize=8.5, loc='left')

ax[1].plot(d/MHz, mn/MHz, color=C[0])
ax[1].axhline(0, color='0.5', lw=0.8)
ax[1].axvspan(-5, u/MHz, color=C[3], alpha=0.15)
ax[1].axvspan(u/MHz, dp/MHz, color=C[1], alpha=0.18)
ax[1].axvspan(dp/MHz, 5, color=C[2], alpha=0.14)
ax[1].axvline(dm/MHz, color=C[1], ls='--', lw=1.0)
ax[1].text(dm/MHz, -5.6, r'$\delta_-$', color=C[1], fontsize=8, ha='center')
ax[1].text(-2.4, 0.9, 'two-photon resonance\n$\\omega_-<0$: not a minimum,\nno symmetry breaking',
           fontsize=6.2, ha='center', color='0.25')
ax[1].text(2.40, -3.6, 'superradiant\nwindow', fontsize=6.4, ha='center', color=C[1])
ax[1].text(3.9, 0.75, 'normal\n(dark)', fontsize=6.4, ha='center', color=C[2])
ax[1].set_xlabel(r'$\delta/2\pi$  (MHz)')
ax[1].set_ylabel(r'min eig $V$  (MHz)')
ax[1].set_ylim(-6.5, 1.6)
ax[1].set_title('(b) only the upper edge is critical', fontsize=8.5, loc='left')

fig.savefig('figs/fig7_roots.pdf'); fig.savefig('figs/fig7_roots.png', dpi=185)
plt.close(fig)
print(f"delta_- = {dm/MHz:+.4f} MHz, delta_+ = {dp/MHz:+.4f} MHz, u = {u/MHz:.4f} MHz")
print("fig7 done")
