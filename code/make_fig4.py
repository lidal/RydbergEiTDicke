"""Fig 4: Dicke-Ising phase diagram, threshold shift and hysteresis."""
import sys, os, json, warnings
warnings.filterwarnings('ignore'); sys.path.insert(0, os.path.dirname(__file__))
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

plt.rcParams.update({'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'figure.dpi': 200,
    'axes.grid': True, 'grid.alpha': 0.25, 'grid.linewidth': 0.5,
    'lines.linewidth': 1.4, 'savefig.bbox': 'tight', 'axes.axisbelow': True})
C = ['#1b4965', '#c1666b', '#4d908e', '#e09f3e']
D = json.load(open('data/ising.json'))

# --- precise threshold-shift check on a fine lambda grid (vectorised) -------
NG = 481
t = np.linspace(0, np.pi, NG); TA, TB = np.meshgrid(t, t, indexing='ij')
SA, SB, CA, CB = np.sin(TA), np.sin(TB), np.cos(TA), np.cos(TB)
NA, NB = 0.5*(1-CA), 0.5*(1-CB)
def alpha(w0, lam, V, w=1.0, z=4):
    al = -(lam/(2*w))*(SA+SB)
    E = w*al**2 - 0.25*w0*(CA+CB) + lam*al*(SA+SB) + 0.5*V*z*NA*NB
    return abs(al[np.unravel_index(np.argmin(E), E.shape)])

fine = np.linspace(0.480, 0.560, 161)
print("Threshold at w0 = +1 on a fine grid (lam_c^(V=0) = 0.5):")
shift = {}
for V in [0.0, 0.5, 2.0, 5.0]:
    a = np.array([alpha(1.0, L, V) for L in fine])
    on = float(fine[int(np.argmax(a > 1e-3))])
    shift[str(V)] = on
    print(f"   V = {V:4.1f}  ->  onset lam = {on:.4f}   (shift {on-0.5:+.4f})")

fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
fig.subplots_adjust(wspace=0.36)

# (a) phase diagram
w0s, lam2 = np.array(D['w0s']), np.array(D['lam2'])
order = {'NORMAL': 0, 'SOLID': 1, 'SRS': 2, 'SR': 3}
M = np.array([[order[p] for p in row] for row in D['labels']]).T
cmap = ListedColormap(['#eef1f4', C[2], C[3], C[1]])
ax[0].pcolormesh(w0s, lam2, M, cmap=cmap, norm=BoundaryNorm([-.5,.5,1.5,2.5,3.5], 4),
                 shading='auto')
ax[0].plot(w0s[w0s > 0], 0.5*np.sqrt(w0s[w0s > 0]), 'k--', lw=1.1)
for lbl, x, y, c in [('NORMAL', 1.1, 0.25, '0.3'), ('SOLID', -1.4, 0.25, 'w'),
                     ('SRS', -1.4, 0.75, '0.15'), ('SR', 0.9, 1.25, 'w')]:
    ax[0].text(x, y, lbl, ha='center', fontsize=7, color=c, weight='bold')
ax[0].set_xlabel(r'$\tilde\omega_0/\omega$'); ax[0].set_ylabel(r'$\lambda/\omega$')
ax[0].set_title(r'(a) $V=2$, $z=4$', fontsize=8.5, loc='left'); ax[0].grid(False)

# (b) threshold at the operating point is untouched by V
for i, V in enumerate(['0.0', '0.5', '2.0']):
    ax[1].plot(D['lams'], D['cuts'][V]['alpha'], color=plt.cm.viridis(i/3.2),
               label=f'$V={V}$')
ax[1].axvline(0.5, color='0.4', ls='--', lw=0.9)
ax[1].set_xlim(0.3, 1.1); ax[1].set_xlabel(r'$\lambda/\omega$')
ax[1].set_ylabel(r'$|\alpha|$'); ax[1].legend(frameon=False)
ax[1].set_title(r'(b) cut at $\tilde\omega_0=+\omega$', fontsize=8.5, loc='left')

# (c) hysteresis
h = D['hysteresis'] if 'lam' in D['hysteresis'] else list(D['hysteresis'].values())[0]
ax[2].plot(h['lam'], h['up'], color=C[0], label='up sweep')
ax[2].plot(h['lam'], h['down'], color=C[1], ls='--', label='down sweep')
ax[2].set_xlabel(r'$\lambda/\omega$'); ax[2].set_ylabel(r'$|\alpha|$')
ax[2].legend(frameon=False)
ax[2].set_title(r'(c) $\tilde\omega_0=-0.5\,\omega$, $V=2$', fontsize=8.5, loc='left')
fig.savefig('figs/fig4_ising.pdf'); plt.close(fig)
print(f"\nhysteresis gap = {h.get('gap', h.get('max_gap')):.4f}")
json.dump(shift, open('data/ising_shift.json', 'w'), indent=1)
print('fig4 done')
