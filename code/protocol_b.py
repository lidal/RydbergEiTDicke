"""
Protocol B read-out: what you actually see on the detector.

The observable is NOT EIT transparency. It is the light *emitted* by the cavity:
below threshold the cavity is dark, above it the superradiant order parameter
switches on. Scanning the two-photon detuning delta produces a sharp EDGE whose
position is the measurand:

    w0t = delta + Omega_mu/2 ,   lasing/superradiant while  w0t < 4 lambda^2/w
    =>  delta_th = 4 lambda^2/w - Omega_mu/2 ,   d(delta_th)/dE = -d/2hbar .

For contrast we also compute the conventional 4-level ladder EIT-AT transmission
on the same axis, where the same field shows up as a doublet splitting.
"""
import numpy as np

hbar, e, a0 = 1.054571817e-34, 1.602176634e-19, 5.29177210903e-11
K = 1774.82 * e * a0 / hbar          # rad/s per (V/m)
MHz = 2 * np.pi * 1e6


def Omega_mu(E_mVcm):
    """MW Rabi frequency (rad/s) for a field in mV/cm."""
    return K * (E_mVcm * 1e-3 * 1e2)


# ------------------------------------------------ Route II order parameter
def nbar_mf(w, w0t, lam):
    """Mean-field photons per atom of the Dicke model (0 in the normal phase)."""
    w0t = np.asarray(w0t, dtype=float)
    mu = w * w0t / (4 * lam ** 2)
    out = np.where((mu < 1) & (w0t > 0), (lam ** 2 / w ** 2) * (1 - np.clip(mu, 0, 1) ** 2), 0.0)
    return out


def delta_threshold(w, lam, E_mVcm):
    """Lock point of Protocol B, in rad/s."""
    return 4 * lam ** 2 / w - Omega_mu(E_mVcm) / 2


# ------------------------------------------- conventional 4-level ladder EIT-AT
def eit_at_transmission(d2, Om_c, Om_mu, g21, g31, g41, Om_p_small=1.0):
    """
    Weak-probe coherence of the ladder |1>-|2>-|3>-|4>, probe on resonance and the
    two-photon detuning d2 scanned (the standard experimental scan).  Returns the
    normalised absorption Im[rho_21]; the AT doublet is split by Om_mu.
    """
    d2 = np.asarray(d2, dtype=float)
    inner = g41 + 1j * d2
    mid = g31 + 1j * d2 + (Om_mu ** 2 / 4) / inner
    denom = g21 + (Om_c ** 2 / 4) / mid
    return np.imag(1j * Om_p_small / 2 / denom)
