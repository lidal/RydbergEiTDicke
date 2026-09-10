"""
Route II with BOTH microwave-dressed branches retained.

The microwave field dresses |r>,|r'> into |+-> = (|r> +- |r'>)/sqrt2 at
delta +- Omega_mu/2.  The cavity-assisted Raman couples |g> to the *bare* |r>,
and |r> = (|+> + |->)/sqrt2, so it couples to BOTH branches with lambda/sqrt2.

Normal-phase stability (Holstein-Primakoff; exact for the threshold in the
thermodynamic limit).  Potential matrix in quadratures (x_a, x_+, x_-):

    V = [[w, 2g, 2g], [2g, w_+, 0], [2g, 0, w_-]],   g = lambda/sqrt2

det V = 0 gives the threshold

    w w_+ w_- = 4 lambda^2 delta ,     w_pm = delta +- Omega_mu/2 .

Keeping only one branch (as an earlier version of this work did) gives
4 lambda^2 = w w_0, which is wrong in both magnitude and, for the read-out
slope, in sign.
"""
import numpy as np

hbar, e, a0 = 1.054571817e-34, 1.602176634e-19, 5.29177210903e-11
K   = 1774.82 * e * a0 / hbar          # rad/s per (V/m)
MHz = 2 * np.pi * 1e6


def Omega_mu(E_mVcm):
    return K * (E_mVcm * 1e-3 * 1e2)


def potential_matrix(w, wp, wm, lam):
    g = lam / np.sqrt(2.0)
    return np.array([[w, 2*g, 2*g], [2*g, wp, 0.0], [2*g, 0.0, wm]])


def is_normal_stable(w, delta, Om, lam):
    """Normal phase stable iff the potential matrix is positive definite."""
    wp, wm = delta + Om/2, delta - Om/2
    if wm <= 0:
        return False
    return np.min(np.linalg.eigvalsh(potential_matrix(w, wp, wm, lam))) > 0


def lam_c(w, delta, Om):
    """Critical coupling at given delta and microwave field: 4 lam^2 delta = w w_+ w_-."""
    wp, wm = delta + Om/2, delta - Om/2
    if wm <= 0 or delta <= 0:
        return 0.0
    return np.sqrt(w * wp * wm / (4.0 * delta))


def delta_th(w, lam, Om):
    """
    Protocol B lock point: solve w delta^2 - 4 lam^2 delta - w (Om/2)^2 = 0
    for the root with w_- > 0.
    """
    u = Om / 2.0
    a, b, c = w, -4*lam**2, -w*u**2
    disc = b*b - 4*a*c
    if disc < 0:
        return np.nan
    roots = [(-b + np.sqrt(disc)) / (2*a), (-b - np.sqrt(disc)) / (2*a)]
    good = [r for r in roots if r - u > 0]
    return max(good) if good else np.nan


def slope_dth_dOm(w, lam, Om):
    """d(delta_th)/d(Omega_mu).  Asymptotes to +1/2; vanishes as Omega_mu -> 0."""
    d = delta_th(w, lam, Om)
    if not np.isfinite(d):
        return np.nan
    return 0.5 * (w * (Om/2.0)) / (w * d - 2 * lam**2)


def nbar_mf(w, delta, Om, lam):
    """
    Mean-field photons per atom above threshold.  With two branches the
    saturation is shared; we use the standard Dicke form driven by the
    distance of the stability determinant from zero, normalised so that the
    single-branch limit is recovered when one branch dominates.
    """
    lc = lam_c(w, delta, Om)
    if lc <= 0 or lam <= lc:
        return 0.0
    mu = (lc / lam) ** 2
    return (lam ** 2 / w ** 2) * (1 - mu ** 2)
