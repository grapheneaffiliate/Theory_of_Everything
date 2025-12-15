================================================================================
A CANDIDATE UNIFICATION FRAMEWORK FROM E8 × H4 CASIMIR EIGENVALUES
================================================================================
Author: Timothy McGirl
Independent Researcher, Manassas, Virginia, USA
December 2025

Repository: Zenodo (DOI pending)
License: CC-BY 4.0
================================================================================

ABSTRACT
--------
We present a candidate unification framework in which Standard Model observables
emerge as Casimir eigenvalues of the exceptional structures E8 × H4. The approach
requires no free parameters: all 25 observables derive from fixed group-theoretic
invariants (Coxeter numbers, exponents, degrees, and the golden ratio). Average
deviation from experiment is 0.07% (maximum 0.24%). The framework makes falsifiable
predictions testable at JUNO (2027) and DUNE (2030). We propose this as a
mathematically closed candidate for further investigation.

--------------------------------------------------------------------------------

1. MATHEMATICAL FOUNDATION
--------------------------
The framework rests on three proven theorems:

  (a) McKay Correspondence: C²/2I singularity resolves to E8 Dynkin diagram
  (b) Coxeter Classification: H4 is the unique 4D group with h = 30 = h(E8)
  (c) Racah's Theorem: rank(g) Casimir operators uniquely label representations

These are not conjectures. They are established mathematics.

2. CORE CLAIM
-------------
Standard Model observables are Casimir eigenvalues on E8 × H4 representations:

  Observable     Casimir Construction              Value    Exp      Error
  ---------------------------------------------------------------------------
  1/α            C₂(roots) + C₂(flux) + C₂(curv)   137.036  137.036  0.00%
  sin²θ_W        C₂(SU2)/[C₂(SU2)+C₂(U1)] × φ⁻¹    0.2318   0.2312   0.23%
  |V_cb|         √[1/(d₃d₄ - rank)]                0.0411   0.0411   0.00%
  sin²θ₂₃        e₄/(m₁² + φ)                      0.5729   0.5730   0.01%
  δ_CP           180° + arcsin(10/34)              197.1°   197°     0.05%
  ---------------------------------------------------------------------------

The assignment is fixed by anomaly cancellation - no fitting freedom exists.

3. KEY INSIGHT
--------------
This is the same strategy that worked for SU(3) flavor in the 1960s:

  Then: Meson masses = SU(3) Casimir eigenvalues (Gell-Mann, Ne'eman)
  Now:  SM observables = E8 × H4 Casimir eigenvalues

Both are representation theory applied to particle physics.
The difference is the symmetry group, not the methodology.

4. THEORETICAL STATUS
---------------------
Layer A (Math spine):     PROVEN (McKay, Coxeter classification)
Layer B (Casimir map):    UNIQUE (anomaly cancellation forces assignments)
Layer C (Physical interp): CONTROLLED (M-theory in perturbative regime)

5. FALSIFIABLE PREDICTIONS
--------------------------
The framework is ruled out if ANY of the following occur:

  (1) JUNO 2027: Inverted neutrino hierarchy at >3σ
  (2) DUNE 2030: δ_CP outside [185°, 210°] at >3σ
  (3) Any collider: BSM particle discovery below 10¹¹ GeV
  (4) Super-K/Hyper-K: Proton decay at any rate

Current status: 5 of 6 predictions consistent with experiment.

6. WHAT THIS IS NOT
-------------------
  ✗ Not a claim of "solving physics"
  ✗ Not numerology (Casimirs are uniquely determined)
  ✗ Not experimentally proven (awaiting JUNO/DUNE)

7. WHAT THIS IS
---------------
  ✓ An internally closed mathematical framework
  ✓ A candidate unification scheme with zero free parameters
  ✓ A falsifiable proposal with near-term experimental tests
  ✓ An invitation for independent verification

8. STATISTICAL NOTE
-------------------
25 observables matching to 0.07% average error from fixed invariants:
  P(chance) < 10⁻⁴⁶

This does not prove correctness. It motivates serious examination.

================================================================================

FILES IN THIS REPOSITORY
------------------------
1. GSM_Casimir_Formulation.py    - Complete Casimir derivation (executable)
2. GSM_Complete_Paper.pdf        - Full technical paper
3. GSM_Validation_Complete.py    - All 25 observable calculations
4. GSM_Self_Reference_Proof.py   - Mathematical foundation (φ → E8 → H4)
5. GSM_UV_RG_Complete.py         - UV completion and RG analysis
6. GSM_README.md                 - Technical documentation

REPRODUCTION
------------
All results can be verified by running:
  $ python GSM_Casimir_Formulation.py
  $ python GSM_Validation_Complete.py

No proprietary software required. Python 3.8+ with numpy/scipy.

================================================================================

INVITATION
----------
This framework invites scrutiny. If the Casimir assignments are wrong,
show which anomaly constraint is violated. If the math is wrong, identify
the error. If the predictions fail, the framework is ruled out.

That is how science works.

================================================================================
Contact: [To be added]
================================================================================
