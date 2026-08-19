"""
Core routines for the critical-Dicke Rydberg-EIT microwave electrometer.

Effective Hamiltonian (Raman-engineered Dicke model, rotating frame):

    H = w a^dag a + w0 Jz + (2*lam/sqrt(N)) (a + a^dag) Jx

with j = N/2.  The measurand (microwave Rabi frequency Omega_mu = d E_mu / hbar)
enters through the dressed-polariton splitting

    w0 = delta + Omega_mu / 2 .

Critical coupling:  lam_c = sqrt(w * w0) / 2.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


# ----------------------------------------------------------------------
# operators
# ----------------------------------------------------------------------
def spin_ops(N):
    """Collective spin operators in the maximal (j = N/2) manifold."""
    j = N / 2.0
    m = np.arange(j, -j - 1, -1)          # m = j, j-1, ..., -j
    dim = len(m)
    Jz = sp.diags(m)
    # J_+ |j,m> = sqrt(j(j+1)-m(m+1)) |j,m+1>
    off = np.sqrt(j * (j + 1) - m[1:] * (m[1:] + 1))
    Jp = sp.diags(off, offsets=1, shape=(dim, dim))
    Jm = Jp.T.tocsr()
    Jx = 0.5 * (Jp + Jm)
    return Jz.tocsr(), Jx.tocsr(), Jp.tocsr(), Jm.tocsr(), dim


def boson_ops(nmax):
    n = np.arange(nmax + 1)
    a = sp.diags(np.sqrt(n[1:]), offsets=1)
    return a.tocsr(), sp.diags(n).tocsr(), nmax + 1


def dicke_H(N, nmax, w, w0, lam):
    """Full Dicke Hamiltonian on the (nmax+1) x (N+1) product space."""
    a, num, db = boson_ops(nmax)
    Jz, Jx, _, _, ds = spin_ops(N)
    Ib, Is = sp.identity(db, format='csr'), sp.identity(ds, format='csr')
    H = (w * sp.kron(num, Is)
         + w0 * sp.kron(Ib, Jz)
         + (2.0 * lam / np.sqrt(N)) * sp.kron(a + a.T, Jx))
    return H.tocsr()


def parity_mask(N, nmax):
    """Parity Pi = exp(i pi (n + Jz + N/2)); returns boolean mask of the even sector."""
    n = np.arange(nmax + 1)
    exc = np.arange(N + 1)                      # N/2 - m  =  number of flipped spins
    tot = n[:, None] + exc[None, :]
    return (tot.ravel() % 2) == 0


def ground_state(N, nmax, w, w0, lam, k=1, tol=1e-12):
    """Lowest eigenpair(s) of the Dicke Hamiltonian restricted to the even-parity sector."""
    H = dicke_H(N, nmax, w, w0, lam)
    mask = parity_mask(N, nmax)
    idx = np.where(mask)[0]
    He = H[idx][:, idx]
    d = He.shape[0]
    if d < 400:
        ev, evec = np.linalg.eigh(He.toarray())
        ev, evec = ev[:k], evec[:, :k]
    else:
        ev, evec = spla.eigsh(He, k=k, which='SA', tol=tol)
        o = np.argsort(ev)
        ev, evec = ev[o], evec[:, o]
    full = np.zeros((H.shape[0], k))
    full[idx, :] = evec
    return ev, full, idx


def qfi_w0(N, nmax, w, w0, lam, h=None):
    """
    Quantum Fisher information of the Dicke ground state w.r.t. the atomic
    splitting w0, from the fidelity susceptibility
        F_Q = 8 (1 - |<psi(w0)|psi(w0+h)>|) / h^2 .
    Only the even-parity sector matters: d_w0 H = Jz is parity even.
    """
    if h is None:
        h = 1e-4 * max(w0, 1e-3)
    _, v0, _ = ground_state(N, nmax, w, w0 - h, lam)
    _, v1, _ = ground_state(N, nmax, w, w0 + h, lam)
    ov = abs(float(v0[:, 0] @ v1[:, 0]))
    ov = min(ov, 1.0)
    return 8.0 * (1.0 - ov) / (2 * h) ** 2


def photon_number(N, nmax, w, w0, lam):
    _, v, _ = ground_state(N, nmax, w, w0, lam)
    a, num, db = boson_ops(nmax)
    Jz, Jx, _, _, ds = spin_ops(N)
    Is = sp.identity(ds, format='csr')
    Ib = sp.identity(db, format='csr')
    psi = v[:, 0]
    n = float(psi @ (sp.kron(num, Is) @ psi))
    jz = float(psi @ (sp.kron(Ib, Jz) @ psi))
    return n, jz


def gap(N, nmax, w, w0, lam):
    """Full (parity-resolved) low-lying spectrum: returns E1-E0 over both sectors."""
    H = dicke_H(N, nmax, w, w0, lam)
    d = H.shape[0]
    if d < 400:
        ev = np.linalg.eigvalsh(H.toarray())[:3]
    else:
        ev = np.sort(spla.eigsh(H, k=3, which='SA', tol=1e-11)[0])
    return ev[1] - ev[0], ev


# ----------------------------------------------------------------------
# mean field (thermodynamic limit)
# ----------------------------------------------------------------------
def mf_order_parameter(w, w0, lam):
    """Returns (alpha = <a>/sqrt(N), n_bar = <a^dag a>/N, jz = <Jz>/N)."""
    lam_c = 0.5 * np.sqrt(w * w0)
    if lam <= lam_c:
        return 0.0, 0.0, -0.5
    mu = (lam_c / lam) ** 2
    s = np.sqrt(max(0.0, 1 - mu ** 2))
    alpha = -(lam / w) * s
    return alpha, alpha ** 2, -0.5 * mu


def mf_dicke_ising(w, w0, lam, V, ngrid=200001):
    """
    Mean field of the Dicke-Ising ('superradiant solid'-flavoured) model:
    infinite-range Rydberg repulsion V penalising the Rydberg fraction
    p_r = (1 - cos theta)/2, energy per atom

        E/N = w a^2 - (w0/2) cos t + 2 lam a sin t + (V/2) p_r^2 ,
        a   = -(lam/w) sin t .

    Returns (theta*, alpha*, n_bar*).
    """
    t = np.linspace(0.0, np.pi, ngrid)
    al = -(lam / w) * np.sin(t)
    pr = 0.5 * (1 - np.cos(t))
    E = w * al ** 2 - 0.5 * w0 * np.cos(t) + 2 * lam * al * np.sin(t) + 0.5 * V * pr ** 2
    i = int(np.argmin(E))
    return t[i], al[i], al[i] ** 2


# ----------------------------------------------------------------------
# open system: linear response of the driven normal phase
# ----------------------------------------------------------------------
def open_response(w, w0, lam, kappa, gamma, eta=1.0):
    """
    Normal-phase Holstein-Primakoff bosons a, b with cavity drive eta:

        da/dt = -(i w + kappa) a - i lam (b + b^dag) - i eta
        db/dt = -(i w0 + gamma) b - i lam (a + a^dag)

    Steady state solves M v = -i eta e, with v = (a, a*, b, b*).
    Returns dict with alpha, d alpha / d w0, slowest relaxation rate, and the
    dynamical matrix eigenvalues.
    """
    M = np.array([
        [-(1j * w + kappa), 0, -1j * lam, -1j * lam],
        [0, (1j * w - kappa), 1j * lam, 1j * lam],
        [-1j * lam, -1j * lam, -(1j * w0 + gamma), 0],
        [1j * lam, 1j * lam, 0, (1j * w0 - gamma)],
    ], dtype=complex)
    rhs = np.array([-1j * eta, 1j * eta, 0, 0], dtype=complex)
    v = np.linalg.solve(M, -rhs)

    hh = 1e-6 * max(abs(w0), 1e-9)
    def _a(x):
        Mx = M.copy()
        Mx[2, 2] = -(1j * x + gamma)
        Mx[3, 3] = (1j * x - gamma)
        return np.linalg.solve(Mx, -rhs)[0]
    dadw0 = (_a(w0 + hh) - _a(w0 - hh)) / (2 * hh)

    ev = np.linalg.eigvals(M)
    rate = float(np.min(-ev.real))
    return dict(alpha=v[0], dalpha_dw0=dadw0, rate=rate, eigs=ev)


def mf_dicke_ising_2sub(w, w0, lam, V, z=4, nseed=13):
    """
    Two-sublattice mean field of the Dicke-Ising model on a bipartite lattice
    (the model of Zhang et al., PRL 110, 090402 (2013)):

        H = w a^dag a + w0 sum_i S^z_i + (2 lam/sqrt(N))(a+a^dag) sum_i S^x_i
            + V sum_<ij> n_i n_j ,      n_i = S^z_i + 1/2 .

    With <S^z_{A,B}> = -cos(t)/2, <S^x_{A,B}> = sin(t)/2, alpha = <a>/sqrt(N)
    eliminated by its own stationarity condition:

        E/N = w alpha^2 - (w0/4)(cos tA + cos tB) + lam alpha (sin tA + sin tB)
              + (V z/2) nA nB ,     alpha = -(lam/2w)(sin tA + sin tB).

    Minimised from many seeds (the energy landscape is multi-valleyed once V>0).
    """
    from scipy.optimize import minimize

    def energy(t):
        tA, tB = t
        al = -(lam / (2 * w)) * (np.sin(tA) + np.sin(tB))
        nA, nB = 0.5 * (1 - np.cos(tA)), 0.5 * (1 - np.cos(tB))
        return (w * al ** 2 - 0.25 * w0 * (np.cos(tA) + np.cos(tB))
                + lam * al * (np.sin(tA) + np.sin(tB)) + 0.5 * V * z * nA * nB)

    seeds = np.linspace(0.0, np.pi, nseed)
    best, bt = np.inf, None
    for sA in seeds:
        for sB in seeds:
            r = minimize(energy, [sA, sB], method='L-BFGS-B',
                         bounds=[(0, np.pi), (0, np.pi)], tol=1e-14)
            if r.fun < best - 1e-13:
                best, bt = float(r.fun), r.x
    tA, tB = bt
    al = abs(-(lam / (2 * w)) * (np.sin(tA) + np.sin(tB)))
    nA, nB = 0.5 * (1 - np.cos(tA)), 0.5 * (1 - np.cos(tB))
    cdw = abs(nA - nB)
    TOL = 1e-4
    if al > TOL and cdw > TOL:
        ph = 'SRS'
    elif al > TOL:
        ph = 'SR'
    elif cdw > TOL:
        ph = 'SOLID'
    else:
        ph = 'NORMAL'
    return dict(alpha=float(al), cdw=float(cdw), phase=ph, energy=best,
                nbar=float(al ** 2), nA=float(nA), nB=float(nB))
