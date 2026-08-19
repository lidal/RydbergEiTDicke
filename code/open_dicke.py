"""
Linearised quantum-Langevin treatment of the driven, damped Dicke sensor in the
normal phase (Holstein-Primakoff bosons a = cavity, b = collective atomic mode).

Quadratures v = (x_a, p_a, x_b, p_b),  x = (o+o^dag)/sqrt2,  <x^2>_vac = 1/2.

    H = w a^dag a + w0 b^dag b + lam (a+a^dag)(b+b^dag)
      = w (x_a^2+p_a^2)/2 + w0 (x_b^2+p_b^2)/2 + 2 lam x_a x_b   (+const)

    v_dot = A v + f + B u ,   B = diag(sqrt(2k),sqrt(2k),sqrt(2g),sqrt(2g))
    <u(t) u(t')^T> = (1/2) I delta(t-t')        (vacuum, symmetrised)

Cavity output:  x_out = sqrt(2k) x_a - x_in .
The measurand enters as dH = dw0 * b^dag b, i.e. a force  s = dw0 (0,0,pb_bar,-xb_bar).
"""
import numpy as np


def drift(w, w0, lam, kappa, gamma):
    return np.array([
        [-kappa,  w,      0.0,   0.0],
        [-w,     -kappa, -2*lam, 0.0],
        [0.0,     0.0,   -gamma, w0],
        [-2*lam,  0.0,   -w0,   -gamma],
    ])


def lam_crit(w, w0, kappa, gamma, lo=0.0, hi=None, tol=1e-12):
    """Largest lam with max Re eig(A) < 0 (the damped superradiant threshold)."""
    if hi is None:
        hi = 5 * (0.5 * np.sqrt(w * w0) + kappa + gamma + 1.0)
    f = lambda L: np.max(np.linalg.eigvals(drift(w, w0, L, kappa, gamma)).real)
    if f(hi) < 0:
        return hi
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return lo


def steady_state(w, w0, lam, kappa, gamma, drive):
    """Mean quadratures for a coherent cavity drive f = (0, drive, 0, 0)."""
    A = drift(w, w0, lam, kappa, gamma)
    f = np.array([0.0, drive, 0.0, 0.0])
    return -np.linalg.solve(A, f)


def response_and_noise(w, w0, lam, kappa, gamma, drive, phi, nu=0.0,
                       n_th_b=0.0):
    """
    Returns (|T|, S_out, S_w0) at Fourier frequency nu for homodyne angle phi.
      T     : output-quadrature signal transfer per unit dw0   [1/sqrt(s) per rad/s]
      S_out : two-sided symmetrised output-quadrature noise PSD (vacuum = 1/2)
      S_w0  : input-referred PSD of w0 fluctuations            [(rad/s)^2 / Hz]
    n_th_b adds excess (thermal/dephasing) noise on the atomic channel.
    """
    A = drift(w, w0, lam, kappa, gamma)
    vbar = steady_state(w, w0, lam, kappa, gamma, drive)
    xb, pb = vbar[2], vbar[3]

    Anu = (-1j * nu) * np.eye(4) - A
    Ainv = np.linalg.inv(Anu)

    s = np.array([0.0, 0.0, pb, -xb])              # force per unit dw0
    c = np.array([np.cos(phi), np.sin(phi), 0.0, 0.0])
    L = np.sqrt(2 * kappa) * c

    T = L @ (Ainv @ s)                              # signal transfer

    B = np.diag([np.sqrt(2*kappa), np.sqrt(2*kappa),
                 np.sqrt(2*gamma), np.sqrt(2*gamma)])
    G = L @ (Ainv @ B) - c                          # noise transfer (row, complex)
    Su = 0.5 * np.diag([1.0, 1.0, 1.0 + 2*n_th_b, 1.0 + 2*n_th_b])
    S_out = float(np.real(G @ Su @ G.conj()))

    S_w0 = S_out / abs(T) ** 2 if abs(T) > 0 else np.inf
    return abs(T), S_out, S_w0


def slow_rate(w, w0, lam, kappa, gamma):
    """Slowest relaxation rate = sensor response bandwidth (rad/s)."""
    return float(np.min(-np.linalg.eigvals(drift(w, w0, lam, kappa, gamma)).real))


# ------------------------------------------------------------------ checks
if __name__ == "__main__":
    w = w0 = 1.0
    k, g = 0.05, 0.005

    # 1. empty cavity must return exactly the vacuum level 1/2
    T, So, _ = response_and_noise(w, w0, 0.0, k, g, drive=1.0, phi=0.3, nu=0.17)
    print(f"[check] empty-cavity output noise = {So:.12f}  (vacuum = 0.5)")

    # 2. with atoms but below threshold, output must still be vacuum
    #    (a passive, lossless-in-the-measured-mode linear network is noiseless
    #     only if gamma=0; with gamma>0 it exceeds 1/2)
    T, So, _ = response_and_noise(w, w0, 0.2, k, 0.0, drive=1.0, phi=0.3, nu=0.17)
    print(f"[check] lam=0.2, gamma=0 output noise = {So:.12f}  (vacuum = 0.5)")

    lc = lam_crit(w, w0, k, g)
    print(f"[check] damped threshold lam_c = {lc:.6f}  (undamped {0.5*np.sqrt(w*w0):.6f})")
    print(f"[check] max Re eig at 0.999 lam_c = {np.max(np.linalg.eigvals(drift(w,w0,0.999*lc,k,g)).real):+.3e}")
