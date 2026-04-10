# G-Box / CSTRU-24

**A Golay-Fibered, Reversible Substitution Layer for the 32.C.U.B.I.T. Computational Field**

[![Version](https://img.shields.io/badge/version-4.0-blue.svg)](#)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](#)
[![Status](https://img.shields.io/badge/checks-9%2F9%20passing-brightgreen.svg)](#)
[![License](https://img.shields.io/badge/license-proprietary-lightgrey.svg)](#license)

Reference implementation and formal specification for the **G-Box**: a native, dynamic, Golay-[24,12,8]-fibered reversible substitution primitive embedded in the depth-1 macro-cube of the [32.C.U.B.I.T. v7.7 architecture](https://lumenhelix.com). Integrated into the full **CSTRU-24 Key Encapsulation Mechanism** with 256-bit keys (effective entropy 120–160 bits; see open problems).

---

## What is this?

The G-Box is a **reversible cryptographic substitution primitive** that combines:

- A non-linear rotor `ρ₅(x) = 5x (mod 32)` with exact modular inverse `ρ₁₃`
- A Golay `[24,12,8]` damping filter that enforces structural invariance
- A rigid-mode construction that guarantees bijectivity
- Slice-local operations that commute when disjoint

It serves as the substitution layer for the CSTRU-24 KEM, a post-quantum key encapsulation mechanism operating over the 32-state CORE-32 substrate.

This repository is a **unified specification-and-implementation release**: the paper and the verification script are shipped together as one auditable unit. Every empirical claim in the paper is reproduced by running a single Python file.

---

## Quick start

### Requirements
- Python 3.10+
- NumPy

### Run the full verification suite

```bash
python gbox_cstru24_v4_reference.py
```

Expected output: **9/9 modules PASS** in ~2 seconds. No GPU, no proprietary libraries, no restricted data.

### Run in Google Colab

Upload `gbox_cstru24_v4_reference.py` to a Colab notebook and execute. It runs end-to-end in the free tier.

---

## Verification modules

The reference script contains 9 independently verifiable modules. Each prints pass/fail and measured timings.

| Module | What it verifies | Status |
|--------|------------------|--------|
| **A** | Rotor arithmetic: `ρ₅⁻¹ = ρ₁₃` (exhaustive over Z₃₂); v2.0 bugs (`ρ₂₁`, `ρ₁₇`) refuted | ✅ PASS |
| **B** | `ρ₅` orbit structure: 10 orbits (4 fixed, 2 two-cycles, 2 four-cycles, 2 eight-cycles) | ✅ PASS |
| **C** | Golay `[24,12,8]` code: QR construction, self-duality, weight distribution (759/2576/759) | ✅ PASS |
| **D** | `proj₂₄` definition and Golay membership check; random rate ≈ 2⁻¹² | ✅ PASS |
| **E** | Damping simulation: 0 invariant violations across 10,000 random cascades | ✅ PASS |
| **F** | G-Box round-trip fidelity: 10,000/10,000 trials (100%) | ✅ PASS |
| **G** | Non-linearity of `ρ₅` over F₂⁵: 624/1,024 pairs (60.94%) non-linear | ✅ PASS |
| **H** | Disjoint-slice commutativity: 10,000/10,000; overlapping slices fail 89.74% | ✅ PASS |
| **I** | State encoding: 160 info bits → 32 bytes, 62.5% efficiency | ✅ PASS |

---

## Claim taxonomy

Following the [32.C.U.B.I.T. v7.7 methodology](https://lumenhelix.com), every claim in the paper is tagged by class:

- **(a) Verified here** — Architectural definitions confirmed by the companion script. Includes rotor arithmetic, orbit structure, Golay code parameters, damping probabilities, round-trip fidelity, non-linearity, commutativity, and state encoding.
- **(b) Measured here** — Timing numbers produced by the script under stated conditions.
- **(c) Conjectured here** — KEM correctness under commutative construction, IND-CPA security reduction, DCPP hardness assumption, NIST security level claims, lattice-scheme comparisons. **These require separate formal treatment and are explicitly labelled throughout the paper.**

Class (a) and (b) claims are reproducible by running the script. Class (c) claims are open problems.

---

## Repository contents

```
.
├── README.md                        # This file
├── GBox_CSTRU24_v4.md              # Paper (markdown source)
├── GBox_CSTRU24_v4.docx            # Paper (Word format)
├── gbox_cstru24_v4_reference.py    # Reference implementation + verification suite
└── LICENSE                         # License terms
```

---

## Measured performance

All measurements single-threaded Python 3.10+ / NumPy:

| Operation | Mean (μs) | Trials | Module |
|-----------|-----------|--------|--------|
| `ρ₅` (batch, per face) | 0.018 | 100k | A |
| `ρ₁₃` inverse (batch, per face) | 0.010 | 100k | A |
| Golay `[24,12,8]` encode | 0.249 | 100k | C |
| G-Box 32-face round-trip | 0.533 | 100k | F |

GPU, compiled-language, and FPGA acceleration are future work.

---

## Version history

**v4.0 (March 2026)** — Current. Unified spec-and-implementation release. Added verification script with 9 modules. Added claim taxonomy (a/b/c). Corrected Golay P matrix using the quadratic-residue construction. Fresh measured timings.

**v3.0** — Corrected two critical errors from v2.0: (a) rotor inverse was misidentified as `ρ₂₁` and `ρ₁₇` when the correct value is `ρ₁₃` (since 5 × 13 ≡ 1 mod 32); (b) KEM Decapsulation was logically broken. Redesigned KEM with Diffie–Hellman-style commutative construction using disjoint slices.

**v2.0** — Corrected six issues from v1.0: Grover bound (from 2³² to 2¹²⁸), key entropy specification, added KEM correctness theorem, added `proj₂₄` definition, added hardness caveats, added external references.

**v1.0** — Initial draft.

---

## The v2.0 bugs (for posterity)

Two critical errors in v2.0 were caught during v3.0 review and are now computationally refuted by Module A of the verification script:

1. **Wrong rotor inverse.** v2.0 claimed `ρ₂₁` was the inverse of `ρ₅` and called it the "Cord involution" with the false arithmetic `21² ≡ 1 (mod 32)`. Actual values: `5 × 21 ≡ 9 (mod 32)` and `21² ≡ 25 (mod 32)`. A separate instance in Algorithm 3.3 used `ρ₁₇`, which also fails: `5 × 17 ≡ 21 (mod 32)`. The correct inverse is `ρ₁₃`: `5 × 13 = 65 ≡ 1 (mod 32)`. Verified exhaustively.

2. **Broken KEM Decapsulation.** v2.0's Decaps attempted to invert the ciphertext with `sk = (h₁,…,h₈)`, but the ciphertext was produced by Encaps using different key blocks `K₁,…,K₈`. Since `K ≠ sk`, the inversion doesn't recover anything useful. v3.0 redesigned the KEM with a commutative construction.

These are documented here explicitly because the discipline of catching and publishing corrections is what separates a serious research programme from a confident-sounding vision document.

---

## Open problems

See §7 of the paper for the full list. Highlights:

1. **Key-entropy gap.** The disjoint-slice commutativity requirement constrains effective key entropy to 120–160 bits, short of the 256-bit target. The v2.0 NIST Level 5 claim is suspended.
2. **Complete IND-CCA2 proof** with exact reduction tightness.
3. **C\* code specification.** Concrete generation procedure for the random replacement code.
4. **Constant-time rigid mode** for side-channel resistance.
5. **FPGA/ASIC hardware port.**
6. **Independent cryptanalysis of DCPP** as instantiated in the 32.C.U.B.I.T. architecture.

---

## Citation

If you use this work, please cite:

```bibtex
@techreport{phillips2026gbox,
  author  = {Phillips, Christopher Gordon},
  title   = {The G-Box: A Golay-Fibered, Reversible Substitution Layer
             for the 32.C.U.B.I.T. Computational Field (v4.0)},
  institution = {Lumen Helix Solutions},
  year    = {2026},
  month   = {March},
  url     = {https://lumenhelix.com},
}
```

---

## Related work

This repository is part of a larger research programme at Lumen Helix:

- **32.C.U.B.I.T. v7.7** — The host computational field architecture
- **CORE-32 / R.U.B.I.C.** — The reversible computing substrate
- **C.A.U.L.D.R.O.N.** — The quantum integration module
- **Omni-Directional S-Box** — The HITL governance terminal

See [lumenhelix.com](https://lumenhelix.com) for the full portfolio.

---

## Contact

**Christopher Gordon Phillips**
Independent Researcher — Lumen Helix Solutions
📧 chris@oiq.to
🌐 https://lumenhelix.com

---

## License

© 2026 Christopher Gordon Phillips / LumenHelix Research. All rights reserved.

See `LICENSE` for terms.
