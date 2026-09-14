"""Numerical validation of Step 3: adiabatic elimination of |2>.

Single-excitation sector of the RWA Hamiltonian (Eq. 3 of the note).
Basis (|1,n=1>, |2,0>, |3,0>, |4,0>).  The RWA Hamiltonian conserves
N_exc = a^d a + s22 + s33 + s44, so this 4x4 block is exact.
"""
import numpy as np

def H_exact(wc, D, dl, Dmu, g0, Omc, Ommu, Gam=0.0):
    """Exact 4-level block.  Gam != 0 makes it non-Hermitian (|2> decays)."""
    D2 = D - 0.5j*Gam
    return np.array([
        [wc,  g0,     0.0,     0.0    ],
        [g0,  D2,     Omc/2,   0.0    ],
        [0.0, Omc/2,  dl,      Ommu/2 ],
        [0.0, 0.0,    Ommu/2,  dl+Dmu ]], dtype=complex)

def H_eff(wc, D, dl, Dmu, g0, Omc, Ommu, Gam=0.0):
    """Effective 3-level model after eliminating |2>, basis (|1,1>,|3,0>,|4,0>)."""
    D2 = D - 0.5j*Gam
    lam = g0*Omc/(2*D2)
    return np.array([
        [wc - g0**2/D2, -lam,             0.0    ],
        [-lam,           dl - Omc**2/(4*D2), Ommu/2 ],
        [0.0,            Ommu/2,          dl+Dmu ]], dtype=complex)

def compare(**kw):
    """Return (exact P-space eigenvalues, effective eigenvalues, max |error|)."""
    ee = np.linalg.eigvals(H_exact(**kw))
    ef = np.linalg.eigvals(H_eff(**kw))
    # discard the exact eigenvalue that correlates with |2> (the one near Delta)
    ee = np.array(sorted(ee, key=lambda z: abs(z - kw['D'])))[1:]
    ee = np.array(sorted(ee, key=lambda z: z.real))
    ef = np.array(sorted(ef, key=lambda z: z.real))
    return ee, ef, np.max(np.abs(ee - ef))

if __name__ == "__main__":
    base = dict(wc=0.0, dl=0.30, Dmu=0.0, g0=1.0, Omc=5.0, Ommu=2.0)

    print("Hermitian case, error vs detuning  (all in MHz/2pi)")
    print(f"{'Delta':>8} {'lambda':>9} {'max err':>11} {'err/lambda':>11} "
          f"{'Omc^2/Delta^2':>14}")
    prev = None
    for D in [20.0, 50.0, 100.0, 200.0, 400.0, 800.0]:
        ee, ef, err = compare(D=D, **base)
        lam = base['g0']*base['Omc']/(2*D)
        print(f"{D:8.0f} {lam:9.4f} {err:11.3e} {err/lam:11.3e} "
              f"{(base['Omc']/D)**2:14.3e}")
        if prev: print(f"{'':8} {'':9} scaling exponent: "
                       f"{np.log(prev[1]/err)/np.log(D/prev[0]):.2f}")
        prev = (D, err)

    print("\nEigenvalues at Delta = 100 (MHz/2pi)")
    ee, ef, err = compare(D=100.0, **base)
    for a, b in zip(ee, ef):
        print(f"  exact {a.real:+12.8f}   effective {b.real:+12.8f}   "
              f"diff {abs(a-b):.2e}")

    print("\nVerify the error formula  E_exact - E_eff = -E |c2|^2,")
    print("  i.e. the fractional energy error is the residual |2> population.")
    print(f"  {'Delta':>6} {'E':>11} {'|c2|^2':>10} {'measured':>11} "
          f"{'predicted':>11} {'ratio':>7}")
    for D in [50.0, 100.0, 200.0, 400.0]:
        ee = compare(D=D, **base)[0].real
        w, v = np.linalg.eigh(H_eff(D=D, **base).real)
        Vc = np.array([base['g0'], base['Omc']/2, 0.0])     # couplings to |2>
        for k in range(3):
            c2sq = ((Vc @ v[:, k])/(w[k] - D))**2
            meas, pred = ee[k] - w[k], -w[k]*c2sq
            print(f"  {D:6.0f} {w[k]:+11.6f} {c2sq:10.3e} {meas:+11.3e} "
                  f"{pred:+11.3e} {meas/pred:7.4f}")

    print("\nVerify the inherited decay rates (isolate each channel)")
    Gam = 6.0
    print("  gamma_3 = Omc^2 Gam / 4 Delta^2   (set g0=0: |3> only talks to |2>)")
    for D in [50.0, 100.0, 200.0, 400.0]:
        kw = dict(base); kw.update(g0=0.0, Ommu=0.0)
        ev = np.linalg.eigvals(H_exact(D=D, Gam=Gam, **kw))
        ev = sorted(ev, key=lambda z: -abs(z.imag))[:2]     # the |2>/|3> pair
        r = min(ev, key=lambda z: abs(z.real - base['dl']))
        print(f"    Delta={D:5.0f}  meas={-2*r.imag:.6f}  "
              f"pred={base['Omc']**2*Gam/(4*D**2):.6f}")
    print("  kappa_abs = g0^2 Gam / Delta^2    (set Omc=0: cavity only talks to |2>)")
    for D in [50.0, 100.0, 200.0, 400.0]:
        kw = dict(base); kw.update(Omc=0.0, Ommu=0.0)
        ev = np.linalg.eigvals(H_exact(D=D, Gam=Gam, **kw))
        r = min(ev, key=lambda z: abs(z.real - base['wc']))
        print(f"    Delta={D:5.0f}  meas={-2*r.imag:.6f}  "
              f"pred={base['g0']**2*Gam/D**2:.6f}")
