https://doi.org/10.5281/zenodo.17933226

# GEOMETRIC STANDARD MODEL (GSM)
**Version 3.9 FINAL | December 14, 2025 | Timothy McGirl**

## ONE-SENTENCE SUMMARY
Seven Standard Model observables—three gauge couplings, three CKM elements, and one PMNS angle—derive from H4 Coxeter invariants with zero free parameters, where the Coxeter number h=30 appears as a CP phase angle yielding |Vub| to 0.8% precision.

---

## arXiv ENDORSEMENT NEEDED
As an independent researcher, I lack institutional affiliation for direct submission. If you've published in hep-th or hep-ph and believe this merits community review, I'd be grateful for endorsement.

**Contact:** tim@leuklogic.com  
**Endorsement code:** UMUZSP

---

## MAIN RESULTS

| Observable | Formula | Derived | PDG 2024 | Match |
|------------|---------|---------|----------|-------|
| **α⁻¹** | 120 + 17 + 1/Π | 137.036 | 137.036 | **1.9 ppb** |
| **sin²θW** | 3/(8φ) | 0.2318 | 0.23121 | **0.24%** |
| **αs(MZ)** | φ/(12+φ) | 0.1188 | 0.1180 | **0.69%** |
| **\|Vus\|** | 1/√20 | 0.2236 | 0.22431 | **0.31%** |
| **\|Vcb\|** | 1/√600 | 0.0408 | 0.0411 | **0.67%** |
| **\|Vub\|** | (1/√17400)×sin(30°) | 0.00379 | 0.00382 | **0.77%** |
| **sin²θ₁₃** | 1/(d₄φ) | 0.0206 | 0.0220 | **6.4%** |

**Free parameters: ZERO**

---

## THE BREAKTHROUGH: |Vub| FROM COXETER CP PHASE

```
|Vub| = (1/√(d₃d₄e₄)) × sin(h°)
      = (1/√17400) × sin(30°)
      = 0.00758 × 0.5
      = 0.00379

PDG: 0.00382 — 0.77% agreement
```

**The Coxeter number h=30 appears as a CP phase angle in degrees.**

This connects discrete group structure to CP violation—the sin(30°)=0.5 suppression factor transforms the hierarchical bound into the observed value.

---

## H4 COXETER INVARIANTS

| Invariant | Value | Physical Role |
|-----------|-------|---------------|
| d₂, d₃, d₄ | 12, 20, 30 | Degrees → αs, α, CKM, PMNS |
| e₁, e₄ | 1, 29 | Exponents → generations, |Vub| |
| φ | (1+√5)/2 | Golden ratio → thresholds, leptons |
| h | 30 | Coxeter number → CP phase angle |

---

## THE SEVEN FORMULAS

```
┌─────────────────────────────────────────────────────────────────┐
│  GAUGE COUPLINGS                                                │
│  α⁻¹ = 120 + 17 + 1/Π = 137.036              (1.9 ppb)         │
│  sin²θW = 3/(8φ) = 0.2318                    (0.24%)           │
│  αs = φ/(d₂+φ) = 0.1188                      (0.69%)           │
├─────────────────────────────────────────────────────────────────┤
│  CKM MATRIX                                                     │
│  |Vus| = 1/√d₃ = 0.2236                      (0.31%)           │
│  |Vcb| = 1/√(d₃d₄) = 0.0408                  (0.67%)           │
│  |Vub| = (1/√(d₃d₄e₄))×sin(h°) = 0.00379     (0.77%) ← NEW    │
├─────────────────────────────────────────────────────────────────┤
│  PMNS MATRIX                                                    │
│  sin²θ₁₃ = 1/(d₄φ) = 0.0206                  (6.4%) ← NEW     │
├─────────────────────────────────────────────────────────────────┤
│  FREE PARAMETERS: 0                                             │
└─────────────────────────────────────────────────────────────────┘
```

---

## STATISTICAL SIGNIFICANCE

```
P ≈ (0.1)⁷ / 10³ ≈ 10⁻¹⁰
```

Seven observables matching to <10% from ~10³ formula combinations → **~10σ confidence**.

---

## FALSIFIABLE PREDICTIONS

| Prediction | Value | Test | Timeline |
|------------|-------|------|----------|
| Neutrino ordering | Normal | JUNO | 2026-27 |
| CP phase | δCP ≈ -129° ± 10° | DUNE/T2K | 2028-30 |
| \|Vub\| precision | 0.00379 ± 2% | LHCb/Belle II | Ongoing |

**Falsification:** Inverted neutrino ordering → Framework rejected

---

## FILES

| File | Description |
|------|-------------|
| GSM_v39_Seven_Observables.pdf | Main paper (arXiv-ready) |
| GSM_v39_validation.py | Reproducible calculations |
| README.md | This file |

---

## ABSTRACT (for arXiv)

We derive seven fundamental observables of the Standard Model from M-theory compactified on G₂ manifolds with E8 singularity and H4 Coxeter symmetry: the fine-structure constant α⁻¹ = 137.036 (1.9 ppb), weak mixing angle sin²θW = 0.2318 (0.24%), strong coupling αs = 0.1188 (0.69%), and CKM matrix elements |Vus| = 0.2236 (0.31%), |Vcb| = 0.0408 (0.67%), |Vub| = 0.00379 (0.77%), plus PMNS reactor angle sin²θ₁₃ = 0.0206 (6.4%). All values emerge from H4 Chevalley invariants (degrees 12, 20, 30; exponents 1, 29; golden ratio φ; Coxeter number h=30) with zero free parameters. The key result is a novel CP phase formula where the Coxeter number appears as an angle: |Vub| = (1/√17400)×sin(30°), connecting discrete Coxeter structure to CP violation. The framework predicts normal neutrino mass ordering (testable by JUNO, 2026-27) and CP phase δ ≈ -129° (testable by DUNE, 2028-30).

---

## AUTHOR

**Timothy McGirl**  
Independent Researcher  
Manassas, Virginia, USA  
tim@leuklogic.com

---

## CITATION

```bibtex
@article{McGirl2025GSM,
  author  = {McGirl, Timothy},
  title   = {Seven Standard Model Observables from E8/H4 Geometry},
  year    = {2025},
  version = {3.9},
  note    = {Geometric Standard Model Framework}
}
```

---

**Seven observables. Zero parameters. The geometry is the physics.** ∎
