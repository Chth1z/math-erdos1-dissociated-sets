"""Pre-publication spot checks: recompute certificates and small-case claims independently of results/*.json.

Usage: python tools/audit_spotcheck.py
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
import time
from fractions import Fraction
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tensor import K_bound, adjugate_element, conj, cyc_mul, delta, mult_matrix, tensor_certificate  # noqa: E402
from dissociated import is_dissociated, is_dissociated_big, is_Q_dissociated, lift  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ok_all = True
sys.set_int_max_str_digits(0)


def report(name: str, ok: bool, detail: str = "") -> None:
    global ok_all
    ok_all &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))


# 1. Example 20: b=5, Q=16
phi = [17, -8, -8, 0, 0]
psi, N = adjugate_element(conj(phi))
report("Example 20: psi_1 as printed", psi == [62017, 48704, 49792, 53704, 52104], str(psi))
report("Example 20: N_1 = 266321 = 11^2*31*71", N == 266321 == 121 * 31 * 71, str(N))
report("Example 20: psi_1 * conj(phi_1) = N_1", cyc_mul(psi, conj(phi)) == [N, 0, 0, 0, 0])
u1 = [psi[(i + 1) % 5] for i in range(5)]
report("Example 20: u~_1 = (48704, 49792, 53704, 52104, 62017)", u1 == [48704, 49792, 53704, 52104, 62017], str(u1))
report("Example 20: u~_1 is 16-dissociated but not 17-dissociated", is_Q_dissociated(u1, 16) and not is_Q_dissociated(u1, 17))
A15 = lift(u1, 3)
report("Example 20: lift with l=3 (15 elements) is dissociated", len(A15) == 15 and is_dissociated(A15))
report("Example 20: a_1 = 62017, c_1 = 147200 coprime", u1[4] == 62017 and 2 * u1[0] + u1[1] == 147200 and gcd(62017, 147200) == 1)
report("K_1(5) = 424/55", K_bound(5, 1) == Fraction(424, 55), str(K_bound(5, 1)))
report("Delta_{5,1} = 11/16", delta(5, 1) == Fraction(11, 16), str(delta(5, 1)))

# 2. Recompute table rows and compare with certificates.json
rows = json.load(open(ROOT / "results" / "certificates.json"))
targets = {(5, 1), (7, 1), (5, 2), (7, 2), (11, 2), (11, 3), (13, 3)}
for row in rows:
    b, s = int(row["b"]), int(row["s"])
    if (b, s) not in targets:
        continue
    Q, l, d = int(row["Q"]), int(row["l"]), b**s
    t0 = time.time()
    cert = tensor_certificate(b, s, Q)
    K = K_bound(b, s)
    conts = [L.cont for L in cert.levels]
    maxu = cert.u_max()
    f_bound = Fraction(maxu, 2 ** (l * (d - 1)))
    cond_Q = Q >= 8 * K and Q % (2**s) == 0 and 2**l <= Q - K
    stored = float(row["f_bound"])
    sha = hashlib.sha256(str(maxu).encode()).hexdigest()
    report(
        f"row ({b},{s}) Q={Q} l={l}: conts={conts}, hypotheses (Q>=8K, 2^s|Q, 2^l<=Q-K)",
        all(c == 1 for c in conts) and cond_Q,
        f"K={float(K):.2f}, {time.time()-t0:.1f}s",
    )
    report(
        f"row ({b},{s}): exact ratio <= stored upward bound {row['f_bound']}, full sha(max u) matches",
        0 <= Fraction(row["f_bound"]) - f_bound < Fraction(1, 10**6) and sha == row["sha256_max_u"],
    )
    report(
        f"row ({b},{s}): max u <= Q^(d-1) Delta (1+10K/Q)",
        Fraction(maxu) <= Fraction(Q) ** (d - 1) * delta(b, s) * (1 + 10 * K / Q),
        f"max u / (Q^(d-1) Delta) = {float(Fraction(maxu) / (Fraction(Q) ** (d - 1) * delta(b, s))):.6f}",
    )
    if d <= 7:
        # a weaker but feasible sanity check: the lift with 3 doublings (3d elements) must be dissociated
        A = lift(cert.full_u(), 3)
        report(f"row ({b},{s}): lift of u with l=3 ({len(A)} elements) is dissociated (meet in the middle, big ints)", is_dissociated_big(A))

# 2b. Small-case claims of Section 6, items 2-3
for (b, Q, l) in [(5, 16, 3), (5, 24, 4), (5, 40, 5), (7, 24, 3)]:
    cert = tensor_certificate(b, 1, Q)
    u = cert.full_u()
    A = lift(u, l)
    report(f"item 2: b={b}, Q={Q}: level primitive and lift with l={l} ({len(A)} elements) dissociated",
           cert.primitive and 2**l <= Q - K_bound(b, 1) and is_dissociated(A))
for Q in (24,):
    cert = tensor_certificate(5, 1, Q)
    u = cert.full_u()
    t0 = time.time()
    report(f"item 3: b=5, Q={Q}: u~_1 exactly {Q}-dissociated (brute force)",
           is_Q_dissociated(u, Q) and not is_Q_dissociated(u, Q + 1), f"{time.time()-t0:.0f}s")
for Q in (108, 116, 124):
    cert = tensor_certificate(5, 2, Q)
    u = cert.full_u()
    report(f"item 2: (5,2) Q={Q}: conts={[L.cont for L in cert.levels]}, 25-element u dissociated (meet in the middle)",
           is_dissociated_big(u))

# 3. Certification of H(b,s) by counting good t (Section 6, item 6)
def count_good(b: int, s: int, kmax: int, tmax: int) -> list[int]:
    found = [0] * kmax
    for t in range(1, tmax + 1):
        cert = tensor_certificate(b, s, 2 ** (s + 1) * t)
        for k in range(1, kmax + 1):
            if all(L.cont == 1 for L in cert.levels[:k]) and gcd(cert.levels[k - 1].a, cert.levels[k - 1].c) == 1:
                found[k - 1] += 1
    return found


for (b, s, tmax) in [(5, 2, 40), (7, 2, 60), (11, 3, 426)]:
    t0 = time.time()
    degs = []
    deg_theta = 1
    for k in range(1, s):
        degs.append((b - 1) * deg_theta)
        deg_theta += degs[-1]
    need = [2 * D + 1 for D in degs]
    found = count_good(b, s, s - 1, tmax)
    report(
        f"H({b},{s}): good t among 1..{tmax} = {found}, needed > 2D_k = {[n - 1 for n in need]}",
        all(f >= n for f, n in zip(found, need)),
        f"{time.time()-t0:.1f}s",
    )

# 4. Leading-coefficient argument behind deg a_k = deg c_k = (b-1) deg theta_k: adj(conj g) = ((2^b+1)/3) * all-ones
import sympy as sp  # noqa: E402

for b in (5, 7, 11, 13):
    g = [2, -1, -1] + [0] * (b - 3)
    M = sp.Matrix(mult_matrix(conj(g)))
    adjM = M.adjugate()
    expected = sp.ones(b, b) * sp.Rational(2**b + 1, 3)
    report(f"adj(m_{{conj g}}) = ((2^{b}+1)/3) * ones for b={b}", adjM == expected)

print("\nALL PASS" if ok_all else "\nSOME CHECKS FAILED")
