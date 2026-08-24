"""
The U(1) route: driven-dissipative Tavis-Cummings (RWA) with incoherent repump.

In a driven optical cavity there is no chemical potential, so the equilibrium
U(1) polariton condensation of Zhang et al. is realised as a *superradiant
lasing* threshold.  Semiclassical Maxwell-Bloch equations for a homogeneously
broadened collective transition:

    a_dot = -(kappa + i D_c) a + (G^2/g) p           (cavity field)
    p_dot = -(Gamma + i D_a) p + g a s               (collective coherence)
    s_dot = -g_par (s - s0) - 2 g (a* p + c.c.)      (inversion)

with D_c = w_cav - nu, D_a = w0t - nu in a frame rotating at the emission
frequency nu, and G^2 = g^2 N s0 / 2 the collective gain parameter.

The measurand enters through   w0t = delta + Omega_mu/2   as the atom-cavity
detuning  D = w_cav - w0t.
"""
import numpy as np


# ---------------------------------------------------------------- threshold
def threshold_G2(kappa, Gamma, D):
    """
    Lasing threshold.  Setting a_dot = p_dot = 0 for infinitesimal a, p:
        G^2 = (kappa + i D_c)(Gamma + i D_a).
    Im part fixes the emission frequency (cavity pulling), Re part the threshold:
        G^2_th = kappa*Gamma * [1 + D^2/(kappa+Gamma)^2].
    """
    return kappa * Gamma * (1.0 + D ** 2 / (kappa + Gamma) ** 2)


def pulling_frequency(w_cav, w0t, kappa, Gamma):
    """nu = (kappa*w0t + Gamma*w_cav)/(kappa+Gamma):  bad cavity (kappa>>Gamma) locks nu to w0t."""
    return (kappa * w0t + Gamma * w_cav) / (kappa + Gamma)


def pulling_coefficients(kappa, Gamma):
    """(d nu/d w0t, d nu/d w_cav).  The second is the cavity-noise leakage."""
    return kappa / (kappa + Gamma), Gamma / (kappa + Gamma)


def threshold_detuning(G2, kappa, Gamma):
    """|D| at which threshold is crossed, for a given collective gain G^2."""
    r = G2 / (kappa * Gamma)
    return (kappa + Gamma) * np.sqrt(r - 1.0) if r > 1 else np.nan


# ------------------------------------------------- above-threshold steady state
def steady_state(G2, kappa, Gamma, D, g_par, g, s0):
    """
    Above threshold the inversion clamps.  With saturation,
        s   = s0 / (1 + I/I_sat) ,
        G^2(s) = G^2_0 * s/s0 ,
    and lasing requires G^2(s) = G^2_th(D).  Hence the clamped inversion and
    the intracavity intensity follow in closed form:
        s/s0 = G2_th/G2_0 ,     I/I_sat = G2_0/G2_th - 1 .
    Returns (photon number in units of the saturation intensity, s/s0).
    """
    G2th = threshold_G2(kappa, Gamma, D)
    if G2 <= G2th:
        return 0.0, 1.0
    return G2 / G2th - 1.0, G2th / G2


# ------------------------------------------------------ numerical verification
def linear_growth_rate(G2, kappa, Gamma, D):
    """
    Largest real part of the eigenvalues of the linearised (a, p) system in the
    lab frame; > 0 means the non-lasing solution is unstable.  Used to check the
    analytic threshold independently.
        d/dt (a, p) = M (a, p),  M = [[-(kappa - i*0), G2^(1/2)], [G2^(1/2), -(Gamma + i D)]]
    written in the frame of the cavity, so the atom carries the full detuning D.
    """
    M = np.array([[-kappa, np.sqrt(G2)],
                  [np.sqrt(G2), -(Gamma + 1j * D)]], dtype=complex)
    return float(np.max(np.linalg.eigvals(M).real))


# ------------------------------------------------ linewidth and sensitivity
def linewidth(kappa, Gamma, n_ph):
    """
    Modified Schawlow-Townes for a bad-cavity laser: the ordinary ST linewidth
    kappa/(2 n) is reduced by the squared cavity-pulling coefficient, because
    phase noise injected into the cavity is only weakly imprinted on the
    atom-locked emission frequency.
        Delta_omega = (kappa / 2 n) * [Gamma/(kappa+Gamma)]^2      [rad/s, FWHM]
    """
    return (kappa / (2.0 * n_ph)) * (Gamma / (kappa + Gamma)) ** 2


def photons_from_repump(N, w, kappa):
    """Above threshold each repumped atom yields one cavity photon: 2 kappa n = N w / 2."""
    return N * w / (4.0 * kappa)


def sensitivity_pulling(kappa, Gamma, N, w, d_over_hbar):
    """
    Frequency-pulling read-out.  Phase diffusion at rate Delta_omega gives
    sqrt(S_nu) = sqrt(Delta_omega); referred through the pulling coefficient and
    then through w0t = delta + Omega_mu/2, E = hbar Omega_mu / d:

        sqrt(S_E) = (2/d_over_hbar) * sqrt(Delta_omega) / c_at      [V/m /sqrt(Hz)]
    Returns (Delta_omega, n_ph, sqrt(S_E) in V/cm/sqrt(Hz)).
    """
    n = photons_from_repump(N, w, kappa)
    dw = linewidth(kappa, Gamma, n)
    c_at, _ = pulling_coefficients(kappa, Gamma)
    sE = (2.0 / d_over_hbar) * np.sqrt(dw) / c_at
    return dw, n, sE * 1e-2


# --------------------------------------- Rydberg interactions: self-consistency
def rydberg_selfconsistent(G2, kappa, Gamma, D0, chi, Imax=None, ngrid=400001):
    """
    Rydberg interactions shift the dressed transition in proportion to the
    Rydberg population, i.e. to the intracavity intensity:  D(I) = D0 - chi*I.
    Above threshold the gain clamps, giving the self-consistency condition

        I = G2 / G2_th(D0 - chi I) - 1 .

    Positive feedback (and hence bistability) requires the shift to pull the
    system *towards* resonance, i.e. sgn(chi) = sgn(D0).  More than one root is
    the optical bistability that underlies the non-equilibrium transition of
    Ding/Wang -- reached here from the lasing side.
    """
    if Imax is None:
        Imax = max(10.0, 3.0 * abs(D0) / max(abs(chi), 1e-12), 5.0 * G2 / (kappa * Gamma))
    I = np.linspace(0.0, Imax, ngrid)
    f = G2 / threshold_G2(kappa, Gamma, D0 - chi * I) - 1.0 - I
    roots = []
    for i in range(len(I) - 1):
        fa, fb = f[i], f[i + 1]
        if fa == 0.0:
            roots.append(I[i])
        elif fa * fb < 0:
            roots.append(I[i] - fa * (I[i + 1] - I[i]) / (fb - fa))
    # de-duplicate (an exact grid-point zero would otherwise be counted twice)
    out = []
    tol = 3.0 * (I[1] - I[0])
    for r in sorted(roots):
        if r >= 0 and (not out or r - out[-1] > tol):
            out.append(r)
    return out
