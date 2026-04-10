# The G-Box: A Golay-Fibered, Reversible Substitution Layer for the 32.C.U.B.I.T. Computational Field

*Integrated into the Full CSTRU-24 Key Encapsulation Mechanism (256-bit Keys)*

**Christopher Gordon Phillips**

Independent Researcher — Lumen Helix Solutions

https://lumenhelix.com • chris@oiq.to

March 2026 — Version 4.0

---

## ABSTRACT

We present the G-Box, a native, dynamic, Golay-[24,12,8]-fibered reversible substitution primitive embedded in the depth-1 macro-cube of the 32.C.U.B.I.T. v7.7 architecture. Each substitution applies the rotor ρ₅(s) = 5s (mod 32) followed by the sympathetic damping operator Δ, which enforces G₂₄ membership on the 24-face projection of the resulting cube state. The inverse rotor is ρ₁₃ (since 5 × 13 = 65 ≡ 1 mod 32).

Version 4.0 consolidates the formal specification and a complete verified reference implementation into a single self-contained release, following the methodology established in 32.C.U.B.I.T. v7.7 [Phillips 2026a]. The accompanying verification script (gbox_cstru24_v4_reference.py) is executable in any Python 3.10+ environment with NumPy and independently reproduces every empirical claim in this paper. Nine verification modules cover: rotor inverse arithmetic, ρ₅ orbit structure, Golay [24,12,8] code verification, proj₂₄ and membership checking, damping probability and rigid-mode statistics, G-Box round-trip fidelity, non-linearity over F₂⁵, disjoint-slice commutativity, and state encoding. All nine modules pass.

The paper also presents a redesigned CSTRU-24 KEM construction using disjoint-slice commutativity, with an honest assessment of the key-entropy constraint this introduces. The KEM correctness theorem and security reduction are provided as class (c) claims requiring further formal treatment; the underlying algebraic primitives are class (a) verified claims.

**Keywords:** reversible S-box, Golay fiber, 32.C.U.B.I.T., DCPP hardness, CSTRU-24 KEM, post-quantum cryptography, verified reference implementation.

---

## 1. Introduction

The 32.C.U.B.I.T. v7.7 field [Phillips 2026a] provides a depth-3 reversible encoding architecture whose every move is a bijection on Z₃₂. Papers I through IV established the move algebra, benchmarks, Golay damping invariance, and domain validation. The missing piece was a native substitution primitive that injects cryptographic non-linearity without sacrificing reversibility or auditability.

The G-Box fills this gap. It is realized geometrically as a 90° slice twist inside the omni-directional macro-cube. We then embed the G-Box into the full CSTRU-24 KEM, scaling the symmetric key to 256 bits while retaining compact ciphertexts and raycast-governed rollback.

Version 4.0 addresses the gap between specification and verification identified in v7.7's methodology. Previous versions (v2.0, v3.0) described the G-Box and CSTRU-24 KEM algorithms as prose specifications separate from any executable code. Version 4.0 ships a single Python script (gbox_cstru24_v4_reference.py) that serves as the reference implementation, benchmark harness, and verification suite. Every empirical number in this paper is produced by running that script.

### 1.1 Version History

Version 2.0 corrected six issues from v1.0 (Grover bound, key entropy, KEM correctness theorem, proj₂₄ definition, hardness caveats, references). Version 3.0 corrected two critical errors in v2.0: (a) the rotor inverse was misidentified as ρ₂₁ and ρ₁₇ when the correct value is ρ₁₃, and (b) the KEM Decapsulation algorithm was logically broken. Version 4.0 adds: (i) the verification script; (ii) the v7.7-style claim taxonomy; (iii) the corrected P matrix (quadratic-residue construction); (iv) measured timings from the unified codebase.

### 1.2 Claim Taxonomy

Following 32.C.U.B.I.T. v7.7 [Phillips 2026a], statements in this paper fall into three classes:

**(a) Verified here:** Architectural definitions verified by the companion script, including rotor arithmetic, orbit structure, Golay code parameters, weight distributions, damping probabilities, round-trip fidelity, non-linearity, commutativity, and state encoding. Each verified claim cites its script module.

**(b) Measured here:** Timing numbers produced by the script under stated conditions (single-threaded Python 3.10+, NumPy).

**(c) Conjectured here:** KEM correctness under the commutative construction, IND-CPA security reduction, DCPP hardness, NIST security level claims, and comparison with lattice-based schemes. These require separate formal treatment.

Class (a) and (b) claims are reproducible by running gbox_cstru24_v4_reference.py. Class (c) claims are explicitly labelled throughout the paper.

---

## 2. Mathematical Foundations of the G-Box

### 2.1 State Space and the Depth-1 Macro-Cube

The depth-1 macro-cube B₁ is a 32-face structure. We identify its face positions with the index set [32] = {0, 1, …, 31}. A state s ∈ S₁ assigns to each face position i an integer value s(i) ∈ Z₃₂. The full state space is S₁ ≅ Z₃₂³².

Individual face operations such as the rotor ρ₅ act slice-locally: when we say ρ₅ is applied to face i of state s, we mean the map s ↦ s' where s'(i) = 5·s(i) (mod 32) and s'(j) = s(j) for all j ≠ i.

**State encoding convention. [Verified: Module I]** Each face value s(i) ∈ Z₃₂ = {0,…,31} requires 5 bits. For byte-aligned storage, we encode each face value in one byte (the 3 high-order bits are zero-padded). A full depth-1 state occupies 32 bytes = 256 bits of storage, of which 160 bits carry information (storage efficiency: 62.5%). All "32-byte" size claims in this paper refer to this byte-aligned encoding.

### 2.2 Rotor Arithmetic [Verified: Module A]

The forward rotor is ρ₅(x) = 5x (mod 32). Its multiplicative inverse is ρ₁₃, since 5 × 13 = 65 ≡ 1 (mod 32). The round-trip ρ₁₃(ρ₅(x)) = x holds for all x ∈ Z₃₂ (verified exhaustively).

**Errata from previous versions.** Version 2.0 identified the inverse as ρ₂₁ (Definition 2.5, Theorem 3.4) and separately as ρ₁₇ (Algorithm 3.3). Both are incorrect: 5 × 21 ≡ 9 (mod 32); 5 × 17 ≡ 21 (mod 32). Additionally, v2.0 called ρ₂₁ the "Cord involution" with the false claim 21² ≡ 1 (mod 32); in fact, 21² = 441 ≡ 25 (mod 32). The term "Cord involution" is retired. Note also that ρ₁₃ is NOT an involution: 13² = 169 ≡ 9 (mod 32) ≠ 1.

### 2.3 ρ₅ Orbit Structure on Z₃₂ [Verified: Module B]

The permutation ρ₅ on Z₃₂ decomposes into 10 orbits:

- 4 fixed points: {0}, {8}, {16}, {24}
- 2 two-cycles: {4, 20}, {12, 28}
- 2 four-cycles: {2, 10, 18, 26}, {6, 30, 22, 14}
- 2 eight-cycles: {1, 5, 25, 29, 17, 21, 9, 13}, {3, 15, 11, 23, 19, 31, 27, 7}

Total: 4(1) + 2(2) + 2(4) + 2(8) = 4 + 4 + 8 + 16 = 32 elements. All 32 covered.

### 2.4 The Golay Projection proj₂₄ [Verified: Modules C, D]

**Definition 2.1 (Golay index set).** Fix a public index set I₂₄ ⊆ [32] with |I₂₄| = 24. In the reference implementation I₂₄ = {0,…,11, 16,…,27} (first 12 and third block of 12 face positions, skipping cap positions {12,…,15} and {28,…,31}).

**Definition 2.2 (Golay projection).** For a state s ∈ S₁, define:

> proj₂₄(s) = ( s(i) mod 2 : i ∈ I₂₄ ) ∈ F₂²⁴

That is, we read the parity of each of the 24 designated face values to obtain a binary vector of length 24.

**The extended binary Golay code G₂₄** is the [24, 12, 8] linear code over F₂ with generator matrix G = [I₁₂ | P], where P is the 12×12 matrix constructed via the quadratic-residue method (QR mod 11 = {1, 3, 4, 5, 9}):

```
P =
  0 1 1 1 1 1 1 1 1 1 1 1
  1 1 1 0 1 1 1 0 0 0 1 0
  1 0 1 1 0 1 1 1 0 0 0 1
  1 1 0 1 1 0 1 1 1 0 0 0
  1 0 1 0 1 1 0 1 1 1 0 0
  1 0 0 1 0 1 1 0 1 1 1 0
  1 0 0 0 1 0 1 1 0 1 1 1
  1 1 0 0 0 1 0 1 1 0 1 1
  1 1 1 0 0 0 1 0 1 1 0 1
  1 1 1 1 0 0 0 1 0 1 1 0
  1 0 1 1 1 0 0 0 1 0 1 1
  1 1 0 1 1 1 0 0 0 1 0 1
```

**Verified properties [Module C]:**

- P · Pᵀ ≡ I₁₂ (mod 2): confirmed (self-duality).
- G · Gᵀ ≡ 0 (mod 2): confirmed (code is self-orthogonal).
- All 4,096 codewords enumerated exhaustively.
- Weight distribution: A₀ = 1, A₈ = 759, A₁₂ = 2,576, A₁₆ = 759, A₂₄ = 1 — exact match.
- Minimum nonzero weight = 8: confirmed (minimum distance d = 8).
- Minimum Hamming distance between distinct codewords ≥ 8: confirmed (50,000 sampled pairs).

**Remark 2.3 (Why parity?).** Taking parity discards the "level" information of each face (the 4 higher-order bits of s(i)) and retains only the least-significant bit. The G₂₄ membership check is a lightweight binary-code constraint imposed on the parity pattern of the active faces. The full face values carry the key material; the parity projection is the non-linearity switch.

### 2.5 The G-Box: Damped Mode and Rigid Mode

**Definition 2.4 (G-Box, damped mode — state-level).** Let s ∈ S₁ be a depth-1 cube state and let σ be a designated slice (a subset of face positions within I₂₄). Define:

> Gᴰ(s, σ) = s'  where s'(i) = 5·s(i) mod 32 for i ∈ σ, s'(i) = s(i) otherwise,
>
> provided proj₂₄(s') ∈ G₂₄. If proj₂₄(s') ∉ G₂₄, then Gᴰ(s, σ) = s (Δ-damped; state unchanged).

**Definition 2.5 (G-Box, rigid mode).** In rigid mode Gʳ, we precondition on non-damping. On the restricted domain Ωʳ ⊆ S₁ where damping does not occur, Gʳ applies ρ₅ coordinatewise on the designated slice, and its inverse is ρ₁₃ applied coordinatewise:

> Gʳ⁻¹(s, σ)(i) = 13·s(i) mod 32  for i ∈ σ;    s(i)  otherwise.

**G-Box round-trip [Verified: Module F].** Forward ρ₅ followed by inverse ρ₁₃ recovers the original state exactly: 10,000/10,000 trials, 100% fidelity. Slice-local variant (I₂₄ only): 10,000/10,000 trials, 100%. Batch verification (100k × 32 faces): exact match confirmed.

### 2.6 Damping Probability [Verified: Module E]

|G₂₄| = 2¹² = 4,096 codewords out of |F₂²⁴| = 2²⁴ = 16,777,216 binary vectors. The probability that a uniformly random parity vector lies in G₂₄ is exactly 2⁻¹² ≈ 0.0244%. Thus approximately 99.976% of rotor outputs are Δ-damped. Rigid mode requires on average 4,096 candidate draws per step.

**Damping simulation [Verified: Module E].** Starting from a valid Golay codeword, 10,000 random cascades were applied. Result: 5 accepted (0.05%), 9,995 damped (99.95%), zero invariant violations. Observed acceptance rate is consistent with the theoretical 2⁻¹² expectation.

### 2.7 Hardness Caveat and Code Substitution

**Critical.** G₂₄ admits a known polynomial-time maximum-likelihood decoder [Berlekamp 1974; Conway & Sloane 1988]. The NP-hardness argument for general syndrome decoding does NOT apply to G₂₄.

**Proposed fix.** For the cryptographic layer of CSTRU-24, replace G₂₄ with a randomly sampled [24, 12, 8]-like linear code C* without the Golay structural symmetry. The G-Box's geometric and reversibility properties depend only on the [n, k, d] parameters and the rigid-mode construction, not on the specific parity-check matrix. Decoding random binary linear codes of this rate remains NP-hard on average [Berlekamp, McEliece & van Tilborg 1978]. **(Class (c) claim: the average-case hardness of C* instances is an assumption, not a theorem.)**

The G₂₄ variant is retained for the theoretical and architectural exposition. All security claims in Section 5 assume C* is in use.

### 2.8 Non-Linearity [Verified: Module G]

**Theorem 2.9 (Non-linearity of ρ₅ over F₂⁵).** The rotor ρ₅ satisfies ρ₅(x ⊕ y) ≠ ρ₅(x) ⊕ ρ₅(y) for 624 of the 1,024 input pairs (x, y) ∈ Z₃₂ × Z₃₂ (non-linearity fraction: 60.94%).

**Proof.** Exhaustive computation over all 1,024 pairs. Concrete counterexample: ρ₅(3 ⊕ 4) = ρ₅(7) = 3, while ρ₅(3) ⊕ ρ₅(4) = 15 ⊕ 20 = 27 ≠ 3. The non-linearity arises from carry propagation in modular multiplication, which does not commute with bitwise XOR. ■

### 2.9 Disjoint-Slice Commutativity [Verified: Module H]

**Proposition 2.10.** If slices σ_a and σ_b are disjoint (σ_a ∩ σ_b = ∅), then keyed rotor operations on these slices commute: the result is independent of application order. Verified: 10,000/10,000 trials with random states and random keys.

**Counterexample for overlapping slices.** When σ_a ∩ σ_b ≠ ∅, commutativity fails in 89.74% of trials (8,974/10,000). This confirms that disjointness is a necessary condition.

---

## 3. Full CSTRU-24 KEM with Embedded G-Box

### 3.1 KEM Design (Commutative Construction)

**(Class (c): This entire section describes a conjectural cryptographic construction. The underlying algebraic primitives — rotor arithmetic, commutativity, round-trip fidelity — are class (a) verified. The KEM protocol, correctness theorem, and security reduction require separate formal treatment.)**

The v2.0 KEM was broken (Decaps applied inverse with sk but ciphertext was produced with K ≠ sk). The v3.0 redesign adopted a Diffie–Hellman-style commutative structure using disjoint face slices.

**Definition 3.1 (Keyed G-Box chain).** For key vector K = (K₁,…,K₈) and pairwise-disjoint slices σ₁,…,σ₈, define:

> Φ_K(s) = Gʳ(…Gʳ(Gʳ(s, σ₁, K₁), σ₂, K₂)…, σ₈, K₈)

where Gʳ(s, σ_i, K_i) applies ρ₅(s(j) ⊕ K_i[j]) for j ∈ σ_i.

**Algorithms.** KeyGen: sk ← random, pk = Φ_sk(G_base). Encaps: r ← random, C = Φ_r(G_base), K = H(Φ_r(pk)). Decaps: K = H(Φ_sk(C)). Correctness relies on Φ_sk(Φ_r(G_base)) = Φ_r(Φ_sk(G_base)) by disjoint-slice commutativity — **but this requires sk and r to act on independently disjoint slice families** (see §3.2).

### 3.2 Key-Entropy Constraint (Open Problem)

The disjoint-partition requirement for commutativity means the face positions must be split into TWO independent families of slices. This halves the effective key material per family. With 24 Golay-fiber positions split into two families of 8 slices × 3 faces × 5 bits, the effective key entropy per family is 60 bits, for a total of 120 bits. Using all 32 face positions yields up to 160 bits. **Neither achieves the 256-bit target from v2.0. The NIST Level 5 claim from v2.0 is suspended.**

### 3.3 The v2.0 Construction as a Symmetric Cipher

The v2.0 KEM construction, while broken as a KEM, is a valid symmetric keyed-bijection cipher when both parties share the same key K. Encrypt: Φ_K(plaintext_state). Decrypt: Φ_K⁻¹ using ρ₁₃ in reverse order. This cipher inherits the G-Box's reversibility and auditability. Security analysis is outside the scope of this paper.

---

## 4. Omni S-Box Governance and Raycast Audit

Every G-Box substitution appears as a visible 90° twist in the 3×3×3 digital tessellation. The Omni S-Box terminal [Phillips 2026c] intercepts each twist via raycast, freezes the node, and displays the full packet (U, C, V, M⃗, μ⃗) before and after the rotor step and G₂₄ membership check. The terminal offers instant ρ₁₃ rollback for human-in-the-loop audit.

---

## 5. Security Analysis

**(Class (c) throughout. All claims in this section are conjectural and require separate formal treatment.)**

### 5.1 Quantum Security Bound

For effective key entropy λ bits, Grover search cost is O(2^(λ/2)). The achievable λ depends on the slice partition design (Section 3.2): currently 120–160 bits, yielding Grover costs of 2⁶⁰–2⁸⁰.

### 5.2 Classical Security: DCPP Hardness

The classical preimage problem reduces to DCPP (Definition 2.10) for C*. Average-case hardness of decoding random binary linear codes is the foundational assumption of McEliece-family cryptosystems [McEliece 1978; Bernstein & Lange 2017], but is not proven unconditionally.

### 5.3 IND-CPA Security Reduction Sketch

If DCPP for C* is hard on average, then CSTRU-24 is IND-CPA secure (game-hop argument: replace real shared key with random; distinguishing requires solving DCPP). The Fujisaki–Okamoto transform [Fujisaki & Okamoto 1999, 2013] upgrades to IND-CCA2 in the random oracle model. **A complete formal proof with reduction tightness is deferred.**

### 5.4 Comparison Table (Honest)

| Scheme | Key Bits | Grover Cost | NIST Level | CT (bytes) |
|--------|----------|-------------|------------|------------|
| AES-128 | 128 | 2⁶⁴ | Level 1 | N/A |
| AES-256 | 256 | 2¹²⁸ | Level 5 | N/A |
| Kyber-1024 | ~3168 (pub) | ~2¹⁷⁸ | Level 5 | 1568 |
| CSTRU-24 (this work) | 120–160 (effective) | 2⁶⁰–2⁸⁰ | TBD | 32 |

Note: comparing 32-byte CSTRU-24 ciphertexts to 1568-byte Kyber ciphertexts is potentially misleading. Lattice ciphertexts carry structured noise enabling security reductions; CSTRU-24's compact state carries no such noise. The compact size is a consequence of the small state space, which also constrains security.

---

## 6. Performance [Measured: Script Modules A, C, F]

| Operation | Mean (μs) | Trials | Module |
|-----------|-----------|--------|--------|
| ρ₅ (batch, per face) | 0.018 | 100k | A |
| ρ₁₃ inverse (batch, per face) | 0.010 | 100k | A |
| Golay [24,12,8] encode | 0.249 | 100k | C |
| G-Box 32-face round-trip | 0.533 | 100k | F |

---

## 7. Open Problems

1. **Key entropy gap.** Resolve the 120–160 bit effective entropy vs. 256-bit target. May require abandoning disjoint-partition commutativity for a more sophisticated algebraic structure.
2. **Complete IND-CCA2 proof.** Formal proof with exact reduction tightness.
3. **C* code specification.** Concrete generation procedure for the random replacement code, with parameter justification and minimum distance verification.
4. **Constant-time rigid mode.** Side-channel-resistant implementation.
5. **FPGA/ASIC hardware port.**
6. **Ciphertext size vs. security tradeoff analysis.**
7. **Independent cryptanalysis of DCPP** as instantiated in the 32.C.U.B.I.T. architecture.

---

## 8. Conclusion

The G-Box v4.0 closes the gap between specification and verification for the CSTRU-24 substitution layer. The companion script gbox_cstru24_v4_reference.py verifies the rotor arithmetic (ρ₅⁻¹ = ρ₁₃, exhaustively confirmed), the ρ₅ orbit structure (10 orbits, all 32 elements covered), the Golay [24,12,8] code (weight distribution exact match, self-duality confirmed), the damping invariance (zero violations across 10,000 cascades), G-Box round-trip fidelity (100% across 10,000 trials), non-linearity (60.94% of F₂⁵ pairs), and disjoint-slice commutativity (100% across 10,000 trials).

The architecture's verified contribution is the precise combination of: a non-linear rotor with an exact modular inverse, a Golay-code damping filter that enforces structural invariance, a rigid-mode construction that guarantees bijectivity, and slice-local operations that commute when disjoint. The more ambitious KEM construction, security claims, and NIST-level comparisons remain class (c) open problems — but the algebraic core is now fully verified and publicly auditable.

**Reproducibility statement.** The file gbox_cstru24_v4_reference.py constitutes both the reference implementation and the verification suite. It requires Python 3.10+ and NumPy. No proprietary libraries or restricted data are needed. Execution produces the exact numbers reported in this paper. The script is provided as supplementary material and will be hosted at the LumenHelix research archive.

---

## Appendix A: Generator Matrix of the Extended Binary Golay Code G₂₄

The extended binary Golay code G₂₄ is a [24, 12, 8] self-dual linear code over F₂. Its generator matrix in standard form G = [I₁₂ | P] where P is the 12×12 matrix constructed via the quadratic-residue method (QR mod 11 = {1, 3, 4, 5, 9}):

```
P =
  0 1 1 1 1 1 1 1 1 1 1 1
  1 1 1 0 1 1 1 0 0 0 1 0
  1 0 1 1 0 1 1 1 0 0 0 1
  1 1 0 1 1 0 1 1 1 0 0 0
  1 0 1 0 1 1 0 1 1 1 0 0
  1 0 0 1 0 1 1 0 1 1 1 0
  1 0 0 0 1 0 1 1 0 1 1 1
  1 1 0 0 0 1 0 1 1 0 1 1
  1 1 1 0 0 0 1 0 1 1 0 1
  1 1 1 1 0 0 0 1 0 1 1 0
  1 0 1 1 1 0 0 0 1 0 1 1
  1 1 0 1 1 1 0 0 0 1 0 1
```

**Verified properties [Module C]:** P · Pᵀ = I₁₂ (mod 2). G · Gᵀ = 0 (mod 2). Weight distribution: A₀ = 1, A₈ = 759, A₁₂ = 2,576, A₁₆ = 759, A₂₄ = 1. Minimum distance d = 8.

The parity-check matrix is H = [Pᵀ | I₁₂]. The code is self-dual (G₂₄ = G₂₄^⊥).

---

## References

[Alagic et al. 2020] G. Alagic et al. Status Report on the Second Round of the NIST PQC Standardization Process. NIST IR 8309, 2020.

[Bellare & Rogaway 1993] M. Bellare, P. Rogaway. Random Oracles are Practical. CCS 1993.

[Berlekamp 1974] E.R. Berlekamp. Key Papers in the Development of Coding Theory. IEEE Press, 1974.

[Berlekamp, McEliece & van Tilborg 1978] E.R. Berlekamp, R.J. McEliece, H.C.A. van Tilborg. On the Inherent Intractability of Certain Coding Problems. IEEE Trans. Info. Theory, 24(3):384–386, 1978.

[Bernstein & Lange 2017] D.J. Bernstein, T. Lange. Post-quantum cryptography. Nature, 549:188–194, 2017.

[Conway & Sloane 1988] J.H. Conway, N.J.A. Sloane. Sphere Packings, Lattices and Groups. Springer, 1988.

[Fujisaki & Okamoto 1999, 2013] E. Fujisaki, T. Okamoto. Secure Integration of Asymmetric and Symmetric Encryption Schemes. J. Cryptology, 26(1):80–101, 2013.

[Golay 1949] M.J.E. Golay. Notes on Digital Coding. Proc. IRE, 37:657, 1949.

[Grover 1996] L.K. Grover. A Fast Quantum Mechanical Algorithm for Database Search. STOC 1996, pp. 212–219.

[McEliece 1978] R.J. McEliece. A Public-Key Cryptosystem Based on Algebraic Coding Theory. DSN Progress Report 42-44, JPL, 1978.

[Niederreiter 1986] H. Niederreiter. Knapsack-Type Cryptosystems and Algebraic Coding Theory. Problems of Control and Information Theory, 15(2):159–166, 1986.

[Phillips 2026a] C.G. Phillips. 32.C.U.B.I.T. Computational Field v7.7. LumenHelix Research, 2026.

[Phillips 2026b] C.G. Phillips. 32.C.U.B.I.T. v7.7 Reference Implementation. LumenHelix Research, 2026.

[Phillips 2026c] C.G. Phillips. Omni-Directional S-Box: Raycast Governance and HITL Audit. LumenHelix Research, 2026.

[Shor 1994] P.W. Shor. Algorithms for Quantum Computation. FOCS 1994, pp. 124–134.

---

© 2026 Christopher Gordon Phillips / LumenHelix Research.
chris@oiq.to • https://lumenhelix.com
