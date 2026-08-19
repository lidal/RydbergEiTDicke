# Critical Dicke electrometry

Feasibility study for extending intracavity Rydberg-EIT microwave electrometry
(Peng *et al.*, [JOSA B **35**, 2272 (2018)](https://doi.org/10.1364/JOSAB.35.002272))
into the Dicke superradiant phase transition
(cf. Zhang *et al.*, [PRL **110**, 090402 (2013)](https://doi.org/10.1103/PhysRevLett.110.090402)),
with sensing as the target application.

**Deliverable:** [`paper/dicke_rydberg_sensing.pdf`](paper/dicke_rydberg_sensing.pdf) (14 pp.)

## The idea in one equation

The microwave field dresses the Rydberg pair, setting the effective two-level splitting
`ω̃₀ = δ + Ω_μ/2` with `Ω_μ = dE_μ/ħ`. The Dicke critical coupling therefore *is* a function
of the measurand:

```
λ_c = ½ √(ω ω̃₀)
```

Bias just below threshold and a weak microwave field becomes a large optical signal.

## Main findings

| | result |
|---|---|
| Ground-state QFI at `λ_c` | `∝ N^{4/3}` (fitted 1.307, N ≤ 512) — but cancels once critical slowing down is costed |
| Normal-phase QFI | `∝ ε^{-2}` (fitted −2.0008), `ε = 1 − λ²/λ_c²` |
| Transduction gain | `∝ ε^{-1}`, free of bandwidth cost until `ε ≲ 3×10⁻⁴` |
| Bandwidth | pinned at `(κ+γ)/2 = 27.5 kHz`, then `∝ ε` (strict gain–bandwidth product, 0.78%) |
| Intrinsic squeezing | up to −10.3 dB, self-generated |
| **Input-referred floor** | `√S_ω̃₀ = 0.72 √(γ/n_b)` — independent of `ε` |
| Field sensitivity | 0.18 nV cm⁻¹Hz^(−1/2) at `n_b = 10⁴` (quantum-noise-limited) |
| Rydberg interactions | do **not** shift the threshold at the operating point (`ω̃₀ > 0`) |

The boxed result is the honest verdict: **criticality is a transducer, not a new scaling law.**
The novel, defensible contribution is a *null* (threshold-crossing) read-out in which the
measurand appears as the laser detuning that holds the system at threshold — SI-traceable,
with a discriminant not limited by the EIT linewidth.

## Layout

```
paper/   LaTeX source + compiled PDF
code/    dicke_core.py    exact diagonalisation, mean field, Dicke-Ising
         open_dicke.py    linearised quantum Langevin (validated against vacuum)
         run_*.py         analysis runs -> data/*.json
         make_figs.py     figures 1-3
         make_fig4.py     figure 4
figs/    generated PDF figures
data/    numerical results (JSON)
```

## Reproducing

```bash
pip install numpy scipy matplotlib
python3 code/run_scaling.py     # finite-size scaling  (~10 min)
python3 code/run_sensing.py     # gain / bandwidth / noise
python3 code/run_control.py     # the fixed-n_b control experiment
python3 code/run_numbers.py     # physical-unit sensitivities
python3 code/run_ising.py       # Dicke-Ising phase diagram  (~16 min)
python3 code/make_figs.py && python3 code/make_fig4.py
cd paper && pdflatex dicke_rydberg_sensing.tex
```

## Literature review — primary sources read

The survey in §6 rests on full texts, not abstracts. Read in full: Zhang *et al.*
(arXiv:1207.4238), Ding *et al.* (arXiv:2207.11947), Wang *et al.* (arXiv:2502.19761).
Read as published abstracts with stated results: Peng 2018, Yang 2020, Wang 2023 (JOSA B,
paywalled full text). Novelty checked by arXiv full-text search across
Dicke/superradiant × Rydberg/electrometry/metrology/sensing.

Three findings changed the write-up:

1. **Zhang *et al.* is rotating-wave.** Their coupling `(g/√N)Σ(b†aψ + h.c.)` conserves
   excitation number (hence the chemical potential in their phase diagram), so their
   "generalised Dicke model" is Tavis–Cummings-type and their superradiance is a U(1)
   polariton condensation — *not* the Z₂ transition at `λ_c = ½√(ωω̃₀)`. Adding the
   counter-rotating terms is a real extension, not a re-derivation. This strengthens the
   proposal but had to be stated correctly.
2. **The enabling engineering already exists.** arXiv:2503.13949 gives a periodic-driving
   recipe for the anisotropic Dicke model in cavity-coupled Rydberg arrays with the
   counter-rotating/rotating ratio tunable 0→∞; arXiv:2511.22230 maps the phase diagram by
   QMC (and independently confirms that Rydberg blockade suppresses cavity occupation —
   the `N_eff` warning in §5). The two-sublattice mean field here is a consistency check
   against that QMC, and agrees on the order of all four transitions.
3. **The recommended precursor experiment is already done.** arXiv:2502.20792 demonstrates
   a cavity-enhanced Rydberg superheterodyne receiver with 19 dB gain — so the de-risking
   step is off the critical path and the project can go straight to the Dicke step.

Near neighbours to the null read-out (arXiv:2605.08535, Adler injection pulling near
synchronisation; arXiv:2306.12544, sub-to-superradiant threshold for Ramsey readout) mean
the novelty claim is worded narrowly: what is unoccupied is reading the measurand as the
*optical two-photon detuning that holds a light–matter phase transition at threshold*.
