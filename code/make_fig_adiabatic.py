"""Figure for the Step-3 note: how well the elimination works."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adiabatic_check import H_exact, H_eff, compare

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                     "mathtext.fontset": "cm", "font.family": "serif"})
base = dict(wc=0.0, dl=0.30, Dmu=0.0, g0=1.0, Omc=5.0, Ommu=2.0)
fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9))

# (a) error vs Delta
D = np.logspace(np.log10(10), np.log10(1000), 40)
err, pred = [], []
for d in D:
    ee = compare(D=d, **base)[0].real
    w, v = np.linalg.eigh(H_eff(D=d, **base).real)
    k = int(np.argmax(np.abs(ee - w)))                # worst-case eigenvalue
    err.append(abs(ee[k] - w[k]))
    pred.append(abs(w[k])*(base["g0"]**2*v[0, k]**2
                           + base["Omc"]**2/4*v[1, k]**2)/d**2)
err, pred = np.array(err), np.array(pred)
ax[0].loglog(D, err, "o", ms=3.2, color="#1b4965", label=r"measured $\max|\delta E|$")
ax[0].loglog(D, pred, "-", lw=1.4, color="#c1666b",
             label=r"$E\,\langle VV^\dagger\rangle/\Delta^2$  (no fit)")
ax[0].loglog(D, base["g0"]*base["Omc"]/(2*D), "--", lw=1.0, color="0.55",
             label=r"$\lambda=g_0\Omega_c/2\Delta$")
ax[0].set_xlabel(r"$\Delta/2\pi$  (MHz)"); ax[0].set_ylabel(r"energy error  (MHz$/2\pi$)")
ax[0].legend(frameon=False, fontsize=7.4, loc="lower left")
ax[0].set_title(r"(a)  error falls as $\Delta^{-2}$, coupling only as $\Delta^{-1}$",
                fontsize=8.5, loc="left")
ax[0].grid(alpha=0.25, lw=0.5)

# (b) zoom on the avoided crossing: the gap IS the quantity Step 3 predicts
Dfix = 100.0
lam  = base["g0"]*base["Omc"]/(2*Dfix)
# locate the crossing and the true minimum gap numerically -- the light shift
# sits on |3> alone, so it shifts the crossing and slightly reshapes the dressing
def Pspace(d):
    E = np.linalg.eigvals(H_exact(D=Dfix, **{**base, "dl": d})).real
    return np.sort(sorted(E, key=lambda x: abs(x - Dfix))[1:])
scan = np.linspace(0.8, 1.3, 20001)
gaps = np.array([Pspace(d)[1] - Pspace(d)[0] for d in scan])
dc   = scan[np.argmin(gaps)]
gmin = gaps.min()
E0   = 0.5*(Pspace(dc)[0] + Pspace(dc)[1])
print(f"  crossing at delta/2pi = {dc:.5f} MHz, min gap = {1e3*gmin:.2f} kHz, "
      f"prediction 2g_- = sqrt(2)*lambda = {1e3*np.sqrt(2)*lam:.2f} kHz")
dl = np.linspace(dc-0.07, dc+0.07, 241)
Ee = np.array([np.sort(np.real(np.linalg.eigvals(
        H_exact(D=Dfix, **{**base, "dl": d})))) for d in dl])
Ef = np.array([np.sort(np.real(np.linalg.eigvals(
        H_eff(D=Dfix, **{**base, "dl": d})))) for d in dl])
for k in [0, 1]:
    ax[1].plot(dl-dc, 1e3*(Ee[:, k]-E0), "-", lw=2.8, color="#9fb8c8",
               label="exact 4-level" if k == 0 else None)
    ax[1].plot((dl-dc)[::5], 1e3*(Ef[::5, k]-E0), ".", ms=3.0, color="#1b4965",
               label="effective 3-level" if k == 0 else None)
ax[1].annotate("", xy=(0, 1e3*gmin/2), xytext=(0, -1e3*gmin/2),
               arrowprops=dict(arrowstyle="<->", color="#c1666b", lw=1.1))
ax[1].text(0.009, -52, r"$2g_-=\sqrt{2}\lambda$" + f"\n$={1e3*np.sqrt(2)*lam:.1f}$ kHz",
           color="#c1666b", fontsize=7.6, va="center", ha="left")
ax[1].set_xlabel(r"$(\delta-\delta_{\rm cross})/2\pi$  (MHz)")
ax[1].set_ylabel(r"$E/2\pi$  (kHz)")
ax[1].legend(frameon=False, fontsize=7.4, loc="upper left")
ax[1].set_title(r"(b)  the avoided crossing, $\Delta/2\pi=100$ MHz",
                fontsize=8.5, loc="left")
ax[1].grid(alpha=0.25, lw=0.5)

fig.tight_layout()
fig.savefig("../paper/fig_adiabatic.pdf", bbox_inches="tight")
print("wrote ../paper/fig_adiabatic.pdf")
