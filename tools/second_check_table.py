"""Second-pass cross-check of Table 1 and the numerical constants quoted in the manuscript.

Recomputes, through src/tensor.py (Bareiss adjugates, a different code path from the Cayley-Hamilton
generator), the exact rationals max u / 2^{l(d-1)} and max u / 2^{l(d-1)+1} for the rows with s <= 3 or
(13,4), rounds them UP to six decimals with integer arithmetic, and compares with the numbers printed in
paper/ejc/ejc-main.tex.  Also checks the Delta bounds, the finite-scale loss claims and the Siegel constants.

Usage: python tools/second_check_table.py
"""
from __future__ import annotations

import math
import re
import sys
import time
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tensor import K_bound, delta, tensor_certificate  # noqa: E402

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]
ok_all = True


def report(name: str, ok: bool, detail: str = "") -> None:
    global ok_all
    ok_all &= ok
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))


def ceil6(x: Fraction) -> str:
    n = -((-x.numerator * 10**6) // x.denominator)  # ceil(10^6 x)
    return f"{n // 10**6}.{n % 10**6:06d}"


def approx5(x: Fraction) -> str:
    n = round(x * 10**5)
    return f"{n // 10**5}.{n % 10**5:05d}"


tex = (ROOT / "paper" / "ejc" / "ejc-main.tex").read_text(encoding="utf-8")
rows = re.findall(r"^(\d+) & (\d+) & (\d+) & (\d+) & (\d+) & (\d+) & ([\d.]+) & ([\d.]+) & ([\d.]+)\\\\$", tex, re.M)
report("Table 1 parsed from the manuscript: 15 rows", len(rows) == 15, str(len(rows)))

for b, s, d, l, n, Q, f_txt, N_txt, D_txt in rows:
    b, s, d, l, n, Q = map(int, (b, s, d, l, n, Q))
    report(f"row ({b},{s}): d, n consistent", d == b**s and n == l * d)
    D = delta(b, s)
    report(f"row ({b},{s}): Delta column {D_txt} is the 5-decimal approximation", approx5(D) == D_txt, approx5(D))
    if s >= 4 and b == 17:
        print(f"       row ({b},{s}): exact recomputation skipped here (3e6-bit integers); covered by tools/verify_certificates.py")
        continue
    t0 = time.time()
    cert = tensor_certificate(b, s, Q)
    K = K_bound(b, s)
    conts = [L.cont for L in cert.levels]
    maxu = Fraction(cert.u_max())
    f_exact = maxu / 2 ** (l * (d - 1))
    N_exact = f_exact / 2
    report(f"row ({b},{s}): certificate conts={conts}, Q>=8K ({Q >= 8*K}), 2^l<=Q-K ({2**l <= Q - K})",
           all(c == 1 for c in conts) and Q >= 8 * K and 2**l <= Q - K, f"{time.time()-t0:.1f}s")
    report(f"row ({b},{s}): f column {f_txt} == ceil6(exact) {ceil6(f_exact)}", ceil6(f_exact) == f_txt)
    report(f"row ({b},{s}): N/2^n column {N_txt} == ceil6(exact/2) {ceil6(N_exact)}", ceil6(N_exact) == N_txt)
    z = Fraction(d - 1) * (Q - 2**l) / 2**l
    report(f"row ({b},{s}): finite-scale envelope 0<=z<1 and 1/(1-z)<1.006", 0 <= z < 1 and 1 / (1 - z) < Fraction(1006, 1000),
           f"1/(1-z)={float(1/(1-z)):.6f}")
    report(f"row ({b},{s}): shape bound max u <= Q^(d-1) Delta (1+10K/Q)", maxu <= Fraction(Q) ** (d - 1) * D * (1 + 10 * K / Q))

# Delta bounds quoted in the text
D113, D133 = delta(11, 3), delta(13, 3)
report("Delta_{11,3} < 0.31618 and > 0.31617 (so the old bound 0.31617 was false)",
       D113 < Fraction(31618, 10**5) and D113 > Fraction(31617, 10**5), f"{float(D113):.7f}")
report("Delta_{13,3} < 0.30299", D133 < Fraction(30299, 10**5), f"{float(D133):.7f}")
report("Delta_{13,3} < 0.303 (Corollary 31)", D133 < Fraction(303, 1000))

# (7,3) gives 0.2317 and (11,3) is the first pair below 0.22002 (Section 6 reading)
print("       (7,3) max A/2^n column:", [r[7] for r in rows if r[0] == '7' and r[1] == '3'])

# Siegel constants: C_d >= (1-K/Q)^(d-1) / (Delta (1+10K/Q)) for rows (11,3) and (17,4)
for b, s in ((11, 3), (17, 4)):
    row = next(r for r in rows if int(r[0]) == b and int(r[1]) == s)
    Q, d = int(row[5]), b**s
    K = K_bound(b, s)
    lower = (1 - K / Q) ** (d - 1) / (delta(b, s) * (1 + 10 * K / Q))
    target = Fraction(315, 100) if b == 11 else Fraction(484, 100)
    report(f"Siegel bound C_{d} >= {float(target)} (computed {float(lower):.5f})", lower >= target)

# Lemma 32 numerics: K_1 <= 2(b^2-1)/3 and K_s <= 5 K_1 d for the table's (b,s)
for b, s in sorted({(int(r[0]), int(r[1])) for r in rows}):
    K1, Ks = K_bound(b, 1), K_bound(b, s)
    report(f"Lemma 32 bounds at ({b},{s}): K_1 <= 2(b^2-1)/3 and K_s <= 5 K_1 b^s",
           K1 <= Fraction(2 * (b * b - 1), 3) and Ks <= 5 * K1 * b**s, f"K_1={float(K1):.2f}, K_s={float(Ks):.1f}")

print("\nALL PASS" if ok_all else "\nSOME CHECKS FAILED")
