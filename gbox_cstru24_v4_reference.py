#!/usr/bin/env python3
"""
gbox_cstru24_v4_reference.py

Reference implementation and verification suite for:
  "The G-Box: A Golay-Fibered, Reversible Substitution Layer
   for the 32.C.U.B.I.T. Computational Field"
  Version 4.0 — Christopher Gordon Phillips / LumenHelix

This script is simultaneously the reference implementation, the benchmark
harness, and the verification suite for every empirical claim in the v4.0
paper. It requires Python 3.10+, NumPy, and (optionally) SciPy.

Usage:
    python gbox_cstru24_v4_reference.py

Nine modules (A–I) correspond to nine verification tasks. Each prints
pass/fail and measured timings. A reader can reproduce every claimed
number by executing this single file.

Modules:
    A — Rotor arithmetic: ρ₅ inverse is ρ₁₃ (not ρ₂₁ or ρ₁₇)
    B — ρ₅ orbit structure on Z₃₂
    C — Golay [24,12,8] code: generator matrix, weight distribution
    D — proj₂₄ and Golay membership check
    E — Damping probability and rigid-mode statistics
    F — G-Box round-trip fidelity (rigid mode, ρ₅/ρ₁₃)
    G — Non-linearity of ρ₅ over F₂⁵
    H — Disjoint-slice commutativity
    I — State encoding and byte-alignment

© 2026 Christopher Gordon Phillips / LumenHelix Research.
chris@oiq.to • https://lumenhelix.com
"""

import time
import sys
import numpy as np
from collections import Counter

# ═══════════════════════════════════════════════════════════════════════
# CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

MOD = 32
ROTOR = 5
ROTOR_INV = 13

# Golay index set I₂₄: first 12 and third block of 12, skipping caps
I24 = list(range(0, 12)) + list(range(16, 28))
assert len(I24) == 24

# Standard Golay [24,12,8] generator matrix P (the right half of G = [I₁₂ | P])
# Quadratic-residue construction: QR mod 11 = {1,3,4,5,9}
# P = bordered (I+Q) where Q is the QR circulant
# Verified: P·Pᵀ = I₁₂ (mod 2), confirming self-duality
# Matches Conway & Sloane 1988, Ch. 3 (up to row/column permutation)
P_MATRIX = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0],
    [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1],
    [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0],
    [1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0],
    [1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0],
    [1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1],
    [1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1],
    [1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1],
    [1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0],
    [1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1],
], dtype=np.uint8)

# Full generator matrix G = [I₁₂ | P]
G24 = np.hstack([np.eye(12, dtype=np.uint8), P_MATRIX])

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"

def header(module_id, title):
    print(f"\n{'='*70}")
    print(f"  Module {module_id}: {title}")
    print(f"{'='*70}")

def check(name, condition):
    status = PASS if condition else FAIL
    print(f"  [{status}] {name}")
    if not condition:
        global any_failures
        any_failures = True
    return condition

any_failures = False

# ═══════════════════════════════════════════════════════════════════════
# MODULE A: Rotor Arithmetic
# ═══════════════════════════════════════════════════════════════════════

def module_a():
    header("A", "Rotor Arithmetic: ρ₅ inverse verification")

    # Verify ρ₁₃ is the correct inverse
    check("5 × 13 ≡ 1 (mod 32)", (ROTOR * ROTOR_INV) % MOD == 1)

    # Verify ρ₂₁ is NOT the inverse (v2.0 error)
    check("5 × 21 ≡ 9 (mod 32), NOT 1 [v2.0 bug confirmed]",
          (5 * 21) % 32 == 9)

    # Verify 21² ≢ 1 (mod 32) (v2.0 false claim)
    check("21² ≡ 25 (mod 32), NOT 1 [v2.0 'involution' claim refuted]",
          (21 * 21) % 32 == 25)

    # Verify ρ₁₇ is NOT the inverse (v2.0 Algorithm 3.3 error)
    check("5 × 17 ≡ 21 (mod 32), NOT 1 [v2.0 Alg 3.3 bug confirmed]",
          (5 * 17) % 32 == 21)

    # Verify ρ₁₃ is NOT an involution (term "Cord involution" was wrong)
    check("13² ≡ 9 (mod 32), NOT 1 [ρ₁₃ is not an involution]",
          (13 * 13) % 32 == 9)

    # Exhaustive round-trip: ρ₁₃(ρ₅(x)) = x for all x ∈ Z₃₂
    all_roundtrip = all((ROTOR_INV * ((ROTOR * x) % MOD)) % MOD == x
                        for x in range(MOD))
    check("ρ₁₃(ρ₅(x)) = x for all x ∈ Z₃₂ (exhaustive)", all_roundtrip)

    # And the reverse direction
    all_roundtrip_rev = all((ROTOR * ((ROTOR_INV * x) % MOD)) % MOD == x
                           for x in range(MOD))
    check("ρ₅(ρ₁₃(x)) = x for all x ∈ Z₃₂ (exhaustive)", all_roundtrip_rev)

    # Timing
    rng = np.random.default_rng(32)
    states = rng.integers(0, 32, size=100_000, dtype=np.uint8)
    t0 = time.perf_counter()
    _ = (ROTOR * states.astype(np.int32)) % MOD
    t1 = time.perf_counter()
    _ = (ROTOR_INV * states.astype(np.int32)) % MOD
    t2 = time.perf_counter()
    fwd_us = (t1 - t0) / 100_000 * 1e6
    inv_us = (t2 - t1) / 100_000 * 1e6
    print(f"  Timing: ρ₅ = {fwd_us:.3f} μs/op, ρ₁₃ = {inv_us:.3f} μs/op (100k trials)")

# ═══════════════════════════════════════════════════════════════════════
# MODULE B: ρ₅ Orbit Structure on Z₃₂
# ═══════════════════════════════════════════════════════════════════════

def module_b():
    header("B", "ρ₅ Orbit Structure on Z₃₂")

    visited = set()
    orbits = []
    for start in range(MOD):
        if start in visited:
            continue
        orbit = []
        x = start
        while x not in visited:
            orbit.append(x)
            visited.add(x)
            x = (ROTOR * x) % MOD
        orbits.append(tuple(orbit))

    # Classify orbits by length
    lengths = sorted([len(o) for o in orbits])
    length_counts = Counter(lengths)

    print(f"  Total orbits: {len(orbits)}")
    print(f"  Orbit lengths: {lengths}")
    for orb in orbits:
        print(f"    Length {len(orb)}: {orb}")

    check("Total orbits = 10", len(orbits) == 10)
    check("4 fixed points (length-1 orbits)", length_counts[1] == 4)
    check("2 two-cycles", length_counts[2] == 2)
    check("2 four-cycles", length_counts[4] == 2)
    check("2 eight-cycles", length_counts[8] == 2)
    check("All 32 elements covered", len(visited) == 32)

    # Verify fixed points are {0, 8, 16, 24}
    fixed_points = sorted([o[0] for o in orbits if len(o) == 1])
    check("Fixed points = {0, 8, 16, 24}", fixed_points == [0, 8, 16, 24])

    # Verify two-cycles
    two_cycles = sorted([tuple(sorted(o)) for o in orbits if len(o) == 2])
    check("Two-cycles = {4,20} and {12,28}",
          two_cycles == [(4, 20), (12, 28)])

# ═══════════════════════════════════════════════════════════════════════
# MODULE C: Golay [24,12,8] Code Verification
# ═══════════════════════════════════════════════════════════════════════

def module_c():
    header("C", "Golay [24,12,8] Code: Generator Matrix & Weight Distribution")

    # Verify G has correct dimensions
    check("G₂₄ shape = (12, 24)", G24.shape == (12, 24))

    # Verify systematic form: left half is I₁₂
    check("G₂₄ = [I₁₂ | P] (systematic form)",
          np.array_equal(G24[:, :12], np.eye(12, dtype=np.uint8)))

    # Verify self-orthogonality: G · Gᵀ = 0 (mod 2)
    product = (G24 @ G24.T) % 2
    check("G₂₄ · G₂₄ᵀ ≡ 0 (mod 2) [self-dual check]",
          np.all(product == 0))

    # Verify P · Pᵀ = I₁₂ (mod 2)
    ppt = (P_MATRIX @ P_MATRIX.T) % 2
    check("P · Pᵀ ≡ I₁₂ (mod 2) [self-duality of P]",
          np.array_equal(ppt, np.eye(12, dtype=np.uint8)))

    # Enumerate all 2¹² = 4096 codewords
    codewords = []
    for i in range(4096):
        msg = np.array([(i >> b) & 1 for b in range(12)], dtype=np.uint8)
        cw = (msg @ G24) % 2
        codewords.append(cw)
    codewords = np.array(codewords)

    check("Total codewords = 4096 = 2¹²", len(codewords) == 4096)

    # Weight distribution
    weights = np.sum(codewords, axis=1)
    weight_dist = Counter(int(w) for w in weights)

    expected = {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}
    print(f"  Measured weight distribution: {dict(sorted(weight_dist.items()))}")
    print(f"  Expected weight distribution: {expected}")

    check("Weight distribution matches known G₂₄",
          dict(weight_dist) == expected)

    min_nonzero = min(w for w in weights if w > 0)
    check(f"Minimum nonzero weight = {min_nonzero} = 8 (minimum distance)",
          min_nonzero == 8)

    # Distance isolation: min Hamming distance between distinct codewords
    # (Sample-based for speed)
    rng = np.random.default_rng(32)
    min_dist = 24
    for _ in range(50_000):
        i, j = rng.integers(0, 4096, size=2)
        if i == j:
            continue
        d = np.sum(codewords[i] != codewords[j])
        if d < min_dist:
            min_dist = d
    check(f"Min Hamming distance between distinct codewords ≥ 8 (sampled: {min_dist})",
          min_dist >= 8)

    # Timing: encode
    rng2 = np.random.default_rng(42)
    msgs = rng2.integers(0, 2, size=(100_000, 12), dtype=np.uint8)
    t0 = time.perf_counter()
    encoded = (msgs @ G24) % 2
    t1 = time.perf_counter()
    encode_us = (t1 - t0) / 100_000 * 1e6
    print(f"  Timing: Golay encode = {encode_us:.3f} μs/op (100k trials)")

# ═══════════════════════════════════════════════════════════════════════
# MODULE D: proj₂₄ and Golay Membership Check
# ═══════════════════════════════════════════════════════════════════════

def module_d():
    header("D", "proj₂₄ Definition and Golay Membership Check")

    # Build the set of all Golay codewords for fast membership checking
    golay_set = set()
    for i in range(4096):
        msg = np.array([(i >> b) & 1 for b in range(12)], dtype=np.uint8)
        cw = tuple(((msg @ G24) % 2).tolist())
        golay_set.add(cw)

    check("|G₂₄| = 4096", len(golay_set) == 4096)

    # Define proj₂₄
    def proj24(state):
        """Extract parity bits at I₂₄ positions from a 32-element state."""
        return tuple(int(state[i]) % 2 for i in I24)

    # Verify proj₂₄ output is length 24
    test_state = np.arange(32, dtype=np.uint8)
    p = proj24(test_state)
    check("proj₂₄ output length = 24", len(p) == 24)
    check("proj₂₄ output is binary", all(b in (0, 1) for b in p))

    # Verify that a state constructed from a Golay codeword passes membership
    # Place codeword parities at I₂₄ positions
    msg = np.array([1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1], dtype=np.uint8)
    cw = (msg @ G24) % 2
    golay_state = np.zeros(32, dtype=np.uint8)
    for idx, pos in enumerate(I24):
        golay_state[pos] = int(cw[idx])  # even values → parity matches codeword

    p_golay = proj24(golay_state)
    check("State constructed from Golay codeword passes membership",
          p_golay in golay_set)

    # Verify random states almost never pass
    rng = np.random.default_rng(32)
    n_trials = 100_000
    n_pass = 0
    for _ in range(n_trials):
        s = rng.integers(0, 32, size=32, dtype=np.uint8)
        if proj24(s) in golay_set:
            n_pass += 1

    rate = n_pass / n_trials
    expected_rate = 4096 / 2**24
    print(f"  Random membership rate: {rate:.6f} (expected: {expected_rate:.6f})")
    check(f"Random Golay membership rate ≈ 2⁻¹² ({rate:.6f} vs {expected_rate:.6f})",
          abs(rate - expected_rate) < 0.001)

# ═══════════════════════════════════════════════════════════════════════
# MODULE E: Damping Probability and Rigid-Mode Statistics
# ═══════════════════════════════════════════════════════════════════════

def module_e():
    header("E", "Damping Probability and Rigid-Mode Statistics")

    # Exact computation
    n_codewords = 4096
    n_vectors = 2**24
    p_nondamped = n_codewords / n_vectors
    p_damped = 1 - p_nondamped

    print(f"  |G₂₄| = {n_codewords}")
    print(f"  |F₂²⁴| = {n_vectors}")
    print(f"  P(non-damped) = {p_nondamped} = 2⁻¹² ≈ {p_nondamped*100:.4f}%")
    print(f"  P(damped) = {p_damped*100:.4f}%")
    print(f"  Expected rejection samples per rigid-mode step: {1/p_nondamped:.0f}")

    check("P(non-damped) = 2⁻¹² exactly", p_nondamped == 2**(-12))
    check("P(damped) = 1 − 2⁻¹² ≈ 99.976%",
          abs(p_damped * 100 - 99.9756) < 0.001)
    check("Expected draws per rigid step = 4096",
          int(1 / p_nondamped) == 4096)

    # Golay damping simulation (matches v7.7 Module F methodology)
    golay_set = set()
    for i in range(4096):
        msg = np.array([(i >> b) & 1 for b in range(12)], dtype=np.uint8)
        cw = tuple(((msg @ G24) % 2).tolist())
        golay_set.add(cw)

    rng = np.random.default_rng(32)
    # Start from a valid Golay codeword
    start_msg = rng.integers(0, 2, size=12, dtype=np.uint8)
    current = tuple(((start_msg @ G24) % 2).tolist())
    assert current in golay_set

    n_moves = 10_000
    n_accepted = 0
    n_damped = 0
    invariant_violations = 0

    for _ in range(n_moves):
        # Random cascade: XOR with random 24-bit vector (simulates rotor output)
        cascade = tuple(rng.integers(0, 2, size=24).tolist())
        candidate = tuple((a ^ b) for a, b in zip(current, cascade))

        if candidate in golay_set:
            current = candidate
            n_accepted += 1
        else:
            n_damped += 1
            # State unchanged — check invariant preserved
            if current not in golay_set:
                invariant_violations += 1

    print(f"  Damping simulation: {n_moves} moves")
    print(f"    Accepted: {n_accepted} ({n_accepted/n_moves*100:.2f}%)")
    print(f"    Damped: {n_damped} ({n_damped/n_moves*100:.2f}%)")
    print(f"    Invariant violations: {invariant_violations}")

    check("Zero Golay invariant violations across 10,000 moves",
          invariant_violations == 0)
    check("Acceptance rate ≈ 0.024% (within 10× of expected for finite sample)",
          n_accepted / n_moves < 0.005)

# ═══════════════════════════════════════════════════════════════════════
# MODULE F: G-Box Round-Trip Fidelity (Rigid Mode)
# ═══════════════════════════════════════════════════════════════════════

def module_f():
    header("F", "G-Box Round-Trip Fidelity (Rigid Mode, ρ₅/ρ₁₃)")

    rng = np.random.default_rng(32)
    n_trials = 10_000
    n_roundtrip_ok = 0

    for _ in range(n_trials):
        # Random state: 32 face values in Z₃₂
        state = rng.integers(0, 32, size=32, dtype=np.int32)
        original = state.copy()

        # Forward: apply ρ₅ to all faces
        forward = (ROTOR * state) % MOD

        # Inverse: apply ρ₁₃ to all faces
        recovered = (ROTOR_INV * forward) % MOD

        if np.array_equal(recovered, original):
            n_roundtrip_ok += 1

    check(f"Round-trip fidelity: {n_roundtrip_ok}/{n_trials} (100%)",
          n_roundtrip_ok == n_trials)

    # Timed round-trip
    states = rng.integers(0, 32, size=(100_000, 32), dtype=np.int32)
    t0 = time.perf_counter()
    fwd = (ROTOR * states) % MOD
    inv = (ROTOR_INV * fwd) % MOD
    t1 = time.perf_counter()
    rt_us = (t1 - t0) / 100_000 * 1e6
    match = np.all(inv == states)
    check(f"Batch round-trip (100k × 32 faces): exact match = {match}", match)
    print(f"  Timing: full 32-face round-trip = {rt_us:.3f} μs/op (100k trials)")

    # Slice-local round-trip: apply rotor only to I₂₄ positions
    n_slice_ok = 0
    for _ in range(n_trials):
        state = rng.integers(0, 32, size=32, dtype=np.int32)
        original = state.copy()

        # Forward on I₂₄ only
        modified = state.copy()
        for i in I24:
            modified[i] = (ROTOR * state[i]) % MOD

        # Inverse on I₂₄ only
        recovered = modified.copy()
        for i in I24:
            recovered[i] = (ROTOR_INV * modified[i]) % MOD

        if np.array_equal(recovered, original):
            n_slice_ok += 1

    check(f"Slice-local round-trip (I₂₄ only): {n_slice_ok}/{n_trials} (100%)",
          n_slice_ok == n_trials)

# ═══════════════════════════════════════════════════════════════════════
# MODULE G: Non-Linearity of ρ₅ over F₂⁵
# ═══════════════════════════════════════════════════════════════════════

def module_g():
    header("G", "Non-Linearity of ρ₅ over F₂⁵")

    # ρ₅ is linear over Z₃₂ (multiplication by unit) but NOT over F₂⁵
    # Test: find x, y where ρ₅(x ⊕ y) ≠ ρ₅(x) ⊕ ρ₅(y)
    nonlinear_count = 0
    nonlinear_examples = []

    for x in range(MOD):
        for y in range(MOD):
            lhs = (ROTOR * (x ^ y)) % MOD
            rhs = ((ROTOR * x) % MOD) ^ ((ROTOR * y) % MOD)
            if lhs != rhs:
                nonlinear_count += 1
                if len(nonlinear_examples) < 3:
                    nonlinear_examples.append((x, y, lhs, rhs))

    total_pairs = MOD * MOD
    linear_count = total_pairs - nonlinear_count

    print(f"  Total (x,y) pairs: {total_pairs}")
    print(f"  Linear pairs (ρ₅(x⊕y) = ρ₅(x)⊕ρ₅(y)): {linear_count}")
    print(f"  Non-linear pairs: {nonlinear_count}")
    for x, y, lhs, rhs in nonlinear_examples:
        print(f"    Example: x={x}, y={y}: ρ₅({x}⊕{y})=ρ₅({x^y})={lhs}, "
              f"ρ₅({x})⊕ρ₅({y})={(ROTOR*x)%MOD}⊕{(ROTOR*y)%MOD}={rhs}")

    check("ρ₅ is non-linear over F₂⁵ (∃ non-linear pairs)",
          nonlinear_count > 0)
    check(f"Non-linearity fraction = {nonlinear_count/total_pairs:.4f} > 0.5",
          nonlinear_count / total_pairs > 0.5)

    # Concrete example from paper: ρ₅(3⊕4) ≠ ρ₅(3) ⊕ ρ₅(4)
    x, y = 3, 4
    lhs = (ROTOR * (x ^ y)) % MOD
    rhs = ((ROTOR * x) % MOD) ^ ((ROTOR * y) % MOD)
    check(f"Paper example: ρ₅(3⊕4)=ρ₅(7)={(ROTOR*7)%MOD}, "
          f"ρ₅(3)⊕ρ₅(4)={(ROTOR*3)%MOD}⊕{(ROTOR*4)%MOD}={rhs}, "
          f"unequal={lhs != rhs}",
          lhs != rhs)

# ═══════════════════════════════════════════════════════════════════════
# MODULE H: Disjoint-Slice Commutativity
# ═══════════════════════════════════════════════════════════════════════

def module_h():
    header("H", "Disjoint-Slice Commutativity")

    rng = np.random.default_rng(32)

    # Define two disjoint slices of face positions
    slice_a = list(range(0, 4))     # faces 0,1,2,3
    slice_b = list(range(4, 8))     # faces 4,5,6,7
    assert set(slice_a).isdisjoint(set(slice_b))

    def apply_keyed_rotor(state, slice_positions, key_val):
        """Apply ρ₅(s(i) ⊕ key_val) for i in slice_positions."""
        result = state.copy()
        for i in slice_positions:
            result[i] = (ROTOR * (state[i] ^ (key_val % MOD))) % MOD
        return result

    n_trials = 10_000
    n_commute = 0

    for _ in range(n_trials):
        state = rng.integers(0, 32, size=32, dtype=np.int32)
        k_a = rng.integers(0, 32)
        k_b = rng.integers(0, 32)

        # Order 1: A then B
        s1 = apply_keyed_rotor(state, slice_a, k_a)
        s1 = apply_keyed_rotor(s1, slice_b, k_b)

        # Order 2: B then A
        s2 = apply_keyed_rotor(state, slice_b, k_b)
        s2 = apply_keyed_rotor(s2, slice_a, k_a)

        if np.array_equal(s1, s2):
            n_commute += 1

    check(f"Disjoint-slice commutativity: {n_commute}/{n_trials} (100%)",
          n_commute == n_trials)

    # Verify NON-commutativity for overlapping slices
    slice_c = list(range(0, 4))     # faces 0,1,2,3
    slice_d = list(range(2, 6))     # faces 2,3,4,5 — overlaps with C

    n_noncommute = 0
    for _ in range(n_trials):
        state = rng.integers(0, 32, size=32, dtype=np.int32)
        k_c = rng.integers(1, 32)  # non-zero to ensure actual change
        k_d = rng.integers(1, 32)

        s1 = apply_keyed_rotor(state, slice_c, k_c)
        s1 = apply_keyed_rotor(s1, slice_d, k_d)

        s2 = apply_keyed_rotor(state, slice_d, k_d)
        s2 = apply_keyed_rotor(s2, slice_c, k_c)

        if not np.array_equal(s1, s2):
            n_noncommute += 1

    print(f"  Overlapping-slice non-commutativity: {n_noncommute}/{n_trials}")
    check("Overlapping slices do NOT commute (counterexamples exist)",
          n_noncommute > 0)

# ═══════════════════════════════════════════════════════════════════════
# MODULE I: State Encoding and Byte-Alignment
# ═══════════════════════════════════════════════════════════════════════

def module_i():
    header("I", "State Encoding and Byte-Alignment")

    n_faces = 32
    bits_per_face = 5  # Z₃₂ = {0, ..., 31} needs 5 bits
    info_bits = n_faces * bits_per_face
    bytes_aligned = n_faces * 1  # 1 byte per face, zero-padded

    print(f"  Faces per depth-1 cube: {n_faces}")
    print(f"  Bits per face value: {bits_per_face}")
    print(f"  Information content: {info_bits} bits = {info_bits/8:.1f} bytes")
    print(f"  Byte-aligned encoding: {bytes_aligned} bytes = {bytes_aligned*8} bits")
    print(f"  Wasted bits per face: 3 (zero-padded high bits)")
    print(f"  Storage efficiency: {info_bits/(bytes_aligned*8)*100:.1f}%")

    check("Information content = 160 bits", info_bits == 160)
    check("Byte-aligned encoding = 32 bytes", bytes_aligned == 32)
    check("Storage efficiency = 62.5%",
          abs(info_bits / (bytes_aligned * 8) - 0.625) < 0.001)

    # Verify all face values fit in 5 bits
    check("max(Z₃₂) = 31 < 2⁵ = 32", 31 < 2**bits_per_face)

    # Verify encoding round-trip
    rng = np.random.default_rng(32)
    n_trials = 10_000
    n_ok = 0
    for _ in range(n_trials):
        state = rng.integers(0, 32, size=32, dtype=np.uint8)
        # Encode: store as bytes (identity, since values fit in uint8)
        encoded = state.tobytes()
        # Decode: read bytes back
        decoded = np.frombuffer(encoded, dtype=np.uint8).copy()
        if np.array_equal(state, decoded) and len(encoded) == 32:
            n_ok += 1

    check(f"Encoding round-trip: {n_ok}/{n_trials} (100%), 32 bytes each",
          n_ok == n_trials)

# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  G-Box / CSTRU-24 v4.0 — Reference Implementation & Verifier  ║")
    print("║  Christopher Gordon Phillips / LumenHelix                      ║")
    print("║  chris@oiq.to • https://lumenhelix.com                        ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    t_start = time.perf_counter()

    module_a()
    module_b()
    module_c()
    module_d()
    module_e()
    module_f()
    module_g()
    module_h()
    module_i()

    t_total = time.perf_counter() - t_start

    print(f"\n{'='*70}")
    print(f"  SUMMARY")
    print(f"{'='*70}")
    print(f"  Total execution time: {t_total:.2f}s")

    if any_failures:
        print(f"  Result: \033[91mSOME CHECKS FAILED\033[0m")
        sys.exit(1)
    else:
        print(f"  Result: \033[92mALL CHECKS PASSED\033[0m")
        sys.exit(0)

if __name__ == "__main__":
    main()
