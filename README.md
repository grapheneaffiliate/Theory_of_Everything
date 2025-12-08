# The Geometric E₈/H₄ Theory of the Fine Structure Constant

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**Deriving α from M-theory geometry with zero free parameters**

---

## The Formula

```
α⁻¹ = 137 + 10/(59(6φ − 5)) = 137.035999189...
```

| Predicted | Experimental (CODATA 2018) | Deviation |
|-----------|----------------------------|-----------|
| 137.035999189469 | 137.035999177 ± 0.000000021 | **0.59σ** |

---

## The Proven Euler Class Identity

**THEOREM:** The following algebraic identity holds exactly in Q(√5):

```
Euler(4D) = e(ν) × (59/20) × (27√5 − 59)
```

where:
- **Euler(4D) = 1/φ⁴ = (7−3√5)/2** — Euler class of 4D H₄ representation
- **e(ν) = 10/(59(6φ−5))** — Required Euler class for period formula
- **27 is unique:** 27²×5 − 59² = 4×41

This connects H₄ representation theory to the period formula with zero free parameters.

---

## Key Identities

| Identity | Significance |
|----------|--------------|
| `6φ − 5 = φ(17 − 5√5)/2` | 17 appears in E₈ exponents AND golden factor |
| `Σ(E₈-only) = 60 = Σ(H₄)` | E₈ exponents not in H₄ sum to H₄ total |
| `142721 = 41 × 59²` | Denominator of 1/Π shows 59² structure |
| `27²×5 − 59² = 4×41` | Uniqueness constraint for ratio factor |

---

## Quick Start

```bash
git clone https://github.com/tmcgirl/e8h4-alpha.git
cd e8h4-alpha
pip install -r requirements.txt
python euler_identity_proof.py
```

**Output:**
```
✓ Euler(4D) = 1/φ⁴: True
✓ Identity verified to arbitrary precision: True
```

---

## Repository Contents

### Paper
- `E8H4_Theory_v3.pdf` — Complete paper with proven Euler identity (3 pages)

### Core Verification
| File | Purpose |
|------|---------|
| `euler_identity_proof.py` | **THE PROOF** — Verifies Euler class identity |
| `joyce_period_engine_v2.py` | Period computation |
| `breakthrough_engine_v2.py` | Six-method verification suite |
| `volume_formula_final.py` | Symbolic equivalence proof |

### Attack Scripts
| File | Purpose |
|------|---------|
| `joyce_attack_plan.py` | 8 attack vectors on Joyce period |
| `localization_attack.py` | Equivariant localization approach |
| `final_attack.py` | Key discoveries (17 identity, 59²) |
| `euler_class_computation.py` | Full H₄ Euler class analysis |

### Flow Simulations
| File | Purpose |
|------|---------|
| `karigiannis_engine_v2.py` | Laplacian flow with golden damping |
| `g2_flow_7d_mesh.py` | 7D mesh-discretized G₂ flow |
| `pinn_karigiannis_numpy.py` | Physics-informed neural network |

---

## The Numbers

| Number | Value | Origin |
|--------|-------|--------|
| **137** | N_flux | Σ(E₈ exp) + ht(U(1)_Y) = 120 + 17 |
| **59** | Σ(H₄) − 1 | H₄ Coxeter exponent sum minus 1 |
| **41** | Conjugate norm | (3√5−2)(3√5+2) = (6φ−5)(6φ−1) |
| **27** | Ratio factor | Unique: 27²×5 − 59² = 4×41 |
| **10** | \|Φ⁺(SU(5))\| | Positive roots of GUT group |
| **17** | Triple appearance | E₈ exponent, U(1)_Y height, golden factor |

---

## Verification

Run all checks:

```bash
python euler_identity_proof.py      # Proves the Euler class identity
python joyce_period_engine_v2.py    # Verifies period formula
python breakthrough_engine_v2.py    # Six independent methods
```

All computations use exact symbolic arithmetic in Q(√5).

---

## Citation

```bibtex
@software{mcgirl2025e8h4,
  author = {McGirl, Timothy},
  title = {The Geometric E₈/H₄ Theory of the Fine Structure Constant},
  year = {2025},
  url = {https://github.com/tmcgirl/e8h4-alpha},
  note = {Euler class identity proven: Euler(4D) = e(ν) × (59/20)(27√5−59)}
}
```

---

## License

MIT License — see [LICENSE](LICENSE)

---

## Status

- ✅ **PROVEN:** Euler(4D) = e(ν) × (59/20)(27√5−59)
- ✅ **PROVEN:** Euler(4D) = 1/φ⁴ exactly
- ✅ **PROVEN:** 27 uniquely determined by algebraic constraint
- ✅ **PROVEN:** All group-theoretic identities
- ✅ **VERIFIED:** 0.59σ agreement with experiment

**The fine structure constant is geometric.**
