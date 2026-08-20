# Prior-art matrix: who has built which piece

Component breakdown of the literature around **critical Dicke electrometry** — the proposal
to extend intracavity Rydberg-EIT microwave sensing ([Peng 2018](https://doi.org/10.1364/JOSAB.35.002272))
into the superradiant phase transition ([Zhang 2013](https://doi.org/10.1103/PhysRevLett.110.090402)).

Columns are the setup elements that matter for the proposal. The last two are the discriminators.

**Legend** — ✓ present · — absent · **bold** = the cell that decides the strand

| Column | Meaning |
|---|---|
| EIT ladder | Rydberg-EIT ladder `\|g⟩→\|e⟩→\|r⟩→\|r'⟩` (or two-photon Rydberg excitation) |
| Optical cavity | Atoms inside an optical resonator |
| Collective g√N | Operated in the collective strong-coupling regime |
| Counter-rot. | Counter-rotating light–matter terms present (⇒ Z₂ superradiant transition possible) |
| Transition | A phase transition is present, and of what kind |
| Rydberg int. | Rydberg–Rydberg interactions play a role in the physics |
| MW = measurand | The microwave field is the quantity being measured (not merely a drive) |

---

## A. The cavity + Rydberg electrometry lineage

Everything except counter-rotating terms — so no transition. Stops at collective Rabi
splitting read out spectroscopically.

| Paper | EIT ladder | Optical cavity | Collective g√N | Counter-rot. | Transition | Rydberg int. | MW = measurand | Th/Exp | Headline |
|---|---|---|---|---|---|---|---|---|---|
| **Peng 2018** — JOSA B **35**, 2272 *(seed)* | ✓ | ✓ | — *(weak coupling)* | **—** | **—** | — | ✓ | Th | 8× lower detectable field, >20× resolution |
| **Yang 2020** — JOSA B **37**, 1664 | ✓ | ✓ | ✓ *(collective Rabi splitting)* | **—** | **—** | — | ✓ | Th | 7× over plain EIT; 4-peak transmission |
| **Wang 2023** — JOSA B **40**, 2604 | ✓ | ✓ | ✓ *(strong coupling)* | **—** | **—** | — | ✓ | Th | 196.7×/26.2×; floor 396.5 nV/cm |
| **Jing 2020** — Nat. Phys. **16**, 911 | ✓ | — | — | — | — | — | ✓ *(superhet)* | Exp | 55 nV cm⁻¹Hz^−½ |
| **Liu 2025** — Chin. Phys. Lett. **42**, 053201 | ✓ | — *(**MW** cavity)* | — | — | — | — | ✓ | Exp | 18 dB power sensitivity |
| **arXiv:2502.20792** — cavity superhet | ✓ | ✓ | ✓ | **—** | **—** | — | ✓ *(superhet)* | Exp | 19 dB sensitivity gain |

## B. The critical-electrometry strand

Has a transition and uses it for sensing — but the **interaction-driven** one (optical
bistability), not light–matter.

| Paper | EIT ladder | Optical cavity | Collective g√N | Counter-rot. | Transition | Rydberg int. | MW = measurand | Th/Exp | Headline |
|---|---|---|---|---|---|---|---|---|---|
| **Ding 2022** — Nat. Phys. **18**, 1447 | ✓ | — | — | — | ✓ **NEPT / bistability** | ✓ *(drives it)* | ✓ | Exp | FI ×10³ vs independent particles; 49 nV cm⁻¹Hz^−½ |
| **Wang 2025** — arXiv:2502.19761 | ✓ | ✓ | — | — | ✓ **NEPT / bistability** | ✓ *(drives it)* | ✓ | Exp | classical FI ×100 from cavity; 2.6 nV cm⁻¹Hz^−½ |

## C. The Dicke / phase-structure strand

Has the transition, the cavity and the Rydberg interactions — and senses nothing.

| Paper | EIT ladder | Optical cavity | Collective g√N | Counter-rot. | Transition | Rydberg int. | MW = measurand | Th/Exp | Headline |
|---|---|---|---|---|---|---|---|---|---|
| **Dimer 2007** — PRA **75**, 013804 / **Baden 2014** — PRL **113**, 020408 | — | ✓ | ✓ | ✓ *(Raman)* | ✓ **Z₂ Dicke** | — | **—** | Th / Exp | Engineered Dicke QPT via cavity-assisted Raman |
| **Zhang 2013** — PRL **110**, 090402 *(seed)* | two-photon | ✓ | ✓ | **— (RWA)** | ✓ **U(1) polariton + SRS** | ✓ | **—** | Th (QMC) | Superradiant solid; first-order boundaries |
| **arXiv:2503.13949** — anisotropic engineering | — | ✓ | ✓ | ✓ *(tunable 0→∞)* | ✓ | ✓ | **—** | Th | Periodic-driving recipe for anisotropic Dicke |
| **arXiv:2511.22230** — Dicke–Ising QMC | — | ✓ | ✓ | ✓ | ✓ **SR / SRS / Solid-½** | ✓ | **— *(MW = drive)*** | Th (QMC) | Orders of all transitions; blockade suppresses cavity occupation |
| **Garbe 2020/2022**, **Rams 2018**, **Gietka 2021/22**, **Gyhm 2025** | — | ✓ | ✓ | ✓ | ✓ | — | — *(generic parameter)* | Th | Critical metrology scaling — and its no-gos |

## D. What this proposal adds

| | EIT ladder | Optical cavity | Collective g√N | Counter-rot. | Transition | Rydberg int. | MW = measurand | Th/Exp |
|---|---|---|---|---|---|---|---|---|
| **Critical Dicke electrometry** | ✓ | ✓ | ✓ | ✓ | ✓ **Z₂ Dicke** | optional | ✓ **sets λ_c** | Th |

---

## Reading the matrix

Every individual column is occupied, and each strand is exactly one column short:

- **Strand A** lacks counter-rotating terms ⇒ no transition at all.
- **Strand B** has criticality and sensing, but the critical point is set by atomic density
  and interaction strength — it drifts with temperature and stray fields, is hysteretic, and
  cannot be tuned quickly.
- **Strand C** has the genuine Z₂ Dicke transition and does no sensing. The near-miss is
  arXiv:2511.22230: microwave fields are *in* the setup, but as the drive that engineers the
  model, never as the measured quantity.

**The empty cell that defines the project** is the intersection of the last two columns —
a Z₂ Dicke transition whose critical point *is* the measurand:

```
λ_c = ½ √(ω ω̃₀),    ω̃₀ = δ + Ω_μ/2,    Ω_μ = d E_μ / ħ
```

Two cells worth staring at:

- **Zhang 2013 → "— (RWA)"**. Their coupling `(g/√N)Σ(b†aψ + h.c.)` conserves excitation
  number, so their "generalised Dicke model" is Tavis–Cummings-type and their superradiance
  is a U(1) polariton condensation. The seed paper is one column short of the model the
  proposal needs — the extension is real, not a re-derivation.
- **arXiv:2503.13949 → "✓ (tunable 0→∞)"**. The enabling technology, already published, and
  it fills exactly the cell Zhang 2013 leaves empty, for the Rydberg platform specifically.

The same RWA restriction runs through all of strand A, so **counter-rotating terms are the
single ingredient separating this proposal from that entire lineage.**

---

## References

| Key | Reference |
|---|---|
| Peng 2018 | Y. Peng *et al.*, *Cavity-enhanced microwave electric field measurement using Rydberg atoms*, J. Opt. Soc. Am. B **35**, 2272 (2018). [doi](https://doi.org/10.1364/JOSAB.35.002272) |
| Yang 2020 | A. Yang *et al.*, *Enhanced measurement of microwave electric fields with collective Rabi splitting*, J. Opt. Soc. Am. B **37**, 1664 (2020). |
| Wang 2023 | Y. Wang *et al.*, *Rydberg-atom-based measurements of microwave electric fields with cavity QED*, J. Opt. Soc. Am. B **40**, 2604 (2023). |
| Jing 2020 | M. Jing *et al.*, *Atomic superheterodyne receiver based on microwave-dressed Rydberg spectroscopy*, Nat. Phys. **16**, 911 (2020). [arXiv:1902.11063](https://arxiv.org/abs/1902.11063) |
| Liu 2025 | *Cavity-enhanced Rydberg atom microwave receiver*, Chin. Phys. Lett. **42**, 053201 (2025). [arXiv:2404.06915](https://arxiv.org/abs/2404.06915) |
| — | *Cavity-enhanced Rydberg atomic superheterodyne receiver*. [arXiv:2502.20792](https://arxiv.org/abs/2502.20792) |
| Ding 2022 | D.-S. Ding, Z.-K. Liu, B.-S. Shi, G.-C. Guo, K. Mølmer, C. S. Adams, *Enhanced metrology at the critical point of a many-body Rydberg atomic system*, Nat. Phys. **18**, 1447 (2022). [arXiv:2207.11947](https://arxiv.org/abs/2207.11947) |
| Wang 2025 | Q. Wang *et al.*, *High-precision measurement of microwave electric field by cavity-enhanced critical behavior…*, Sci. China Phys. Mech. Astron. (2025). [arXiv:2502.19761](https://arxiv.org/abs/2502.19761) |
| Dimer 2007 | F. Dimer, B. Estienne, A. S. Parkins, H. J. Carmichael, PRA **75**, 013804 (2007). |
| Baden 2014 | M. P. Baden *et al.*, *Realization of the Dicke model using cavity-assisted Raman transitions*, PRL **113**, 020408 (2014); erratum **118**, 199901 (2017). |
| Zhang 2013 | X.-F. Zhang, Q. Sun, Y.-C. Wen, W.-M. Liu, S. Eggert, A.-C. Ji, *Rydberg polaritons in a cavity: a superradiant solid*, PRL **110**, 090402 (2013). [arXiv:1207.4238](https://arxiv.org/abs/1207.4238) |
| — | *Engineering anisotropic Dicke model with dipole-dipole interaction for Rydberg atom arrays in cavity*. [arXiv:2503.13949](https://arxiv.org/abs/2503.13949) |
| — | *Quantum phase transitions of the anisotropic Dicke–Ising model in driven Rydberg arrays*, PRA (2026). [arXiv:2511.22230](https://arxiv.org/abs/2511.22230) |
| Garbe 2020 | L. Garbe, M. Bina, A. Keller, M. G. A. Paris, S. Felicetti, PRL **124**, 120504 (2020). |
| Garbe 2022 | L. Garbe, O. Abah, S. Felicetti, R. Puebla, Quantum Sci. Technol. **7**, 035010 (2022). [arXiv:2110.04144](https://arxiv.org/abs/2110.04144) |
| Rams 2018 | M. M. Rams *et al.*, *At the limits of criticality-based quantum metrology*, PRX **8**, 021022 (2018). [arXiv:1702.05660](https://arxiv.org/abs/1702.05660) |
| Gyhm 2025 | J.-Y. Gyhm, H. Kwon, M.-J. Hwang, *Fundamental scaling limit in critical quantum metrology*. [arXiv:2506.19003](https://arxiv.org/abs/2506.19003) |

*Full discussion in [`dicke_rydberg_sensing.pdf`](dicke_rydberg_sensing.pdf) §6.*
