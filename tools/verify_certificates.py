"""Independently verify every exact witness in results/certificates.json.

No src module or production adjugate implementation is imported. Adjugates are
reconstructed backwards from two entries, with checked integer divisions. A
separate scalar power-sum recurrence computes the determinant. This establishes
the actual adjugate, rather than checking a scale-invariant identity alone.

Run: python tools/verify_certificates.py [--manifest PATH]
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from functools import cache
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time

sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[1]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def hex_integer(text):
    require(isinstance(text, str) and re.fullmatch(r"[1-9a-f][0-9a-f]*", text), "noncanonical positive hex integer")
    return int(text, 16)


def reconstruct_adjugate(b, theta, rho, seeds):
    """M_ij: a on the diagonal, -theta at j=i+1 and i+2 (cyclic)."""
    require(len(seeds) == 2, "exactly two adjugate seeds are required")
    a = 2 * theta + rho
    psi = [seeds[0], seeds[1]] + [0] * (b - 2)
    for i in range(b - 1, 0, -1):
        value, remainder = divmod(theta * (psi[(i + 1) % b] + psi[(i + 2) % b]), a)
        require(remainder == 0, f"non-exact backwards division at coordinate {i}")
        if i == 1:
            require(value == seeds[1], "the cyclic recurrence does not close")
        psi[i] = value
    # r1,r2 roots of theta*x^2+theta*x-a; T_j=theta^j*(r1^j+r2^j).
    T0, T1 = 2, -theta
    for _ in range(2, b + 1):
        T0, T1 = T1, -theta * T1 + a * theta * T0
    determinant = a ** b + T1 - theta ** b
    require(determinant > 0, "nonpositive determinant")
    require(a * psi[0] - theta * (psi[1] + psi[2]) == determinant,
            "the seed vector is not the genuine adjugate column")
    require(all(x > 0 for x in psi), "nonpositive adjugate entry")
    return psi, determinant


@cache
def independent_K(b, s):
    """Rebuild the base Gram inverse directly, then evaluate the stated norms."""
    import sympy as sp
    B = sp.zeros(b, b - 1)
    for i in range(b - 1):
        B[i, i] += 1
        B[(i + 1) % b, i] -= sp.Rational(1, 2)
        B[(i + 2) % b, i] -= sp.Rational(1, 2)
    D = B * (B.T * B).inv(method="DM")
    duals = [[Fraction(int(D[i,j].p), int(D[i,j].q)) for i in range(b)] for j in range(b - 1)]
    K = sum(abs(x) for col in duals for x in col)
    for level in range(2, s + 1):
        dprime = b ** (level - 1)
        total = Fraction(0)
        for col in duals:
            c = (2 * col[0] + col[1]) / 3
            shift = c * Fraction(dprime - 1, dprime)
            total += sum(abs(x - shift) for x in col) + abs(c) * b * Fraction(dprime - 1, dprime)
        K = dprime * total + Fraction(2 * b, 3) * K
    return K


def verify_row(row, directory):
    b, s, Q, l = (row[key] for key in ("b", "s", "Q", "l"))
    require(row["format"] == "tensor-certificate-v2", "unknown manifest format")
    require(b >= 5 and b % 2 == 1 and b % 3 and s >= 1, "parameters outside theorem scope")
    require(Q > 0 and Q % 2 ** s == 0 and l >= 1, "invalid Q or l")
    require(row["d"] == b ** s and row["n"] == l * b ** s, "incorrect dimension or size")
    witness_path = (directory / row["witness"]).resolve()
    require(witness_path.is_relative_to(directory.resolve()), "witness path leaves the results directory")
    compressed = witness_path.read_bytes()
    require(hashlib.sha256(compressed).hexdigest() == row["sha256_witness_gzip"], "compressed witness hash mismatch")
    payload = gzip.decompress(compressed)
    require(hashlib.sha256(payload).hexdigest() == row["sha256_witness_json"], "witness JSON hash mismatch")
    witness = json.loads(payload)
    require(witness["format"] == "tensor-adjugate-seeds-v1", "unknown witness format")
    require([witness[key] for key in ("b", "s", "Q")] == [b,s,Q], "witness parameters differ")
    require(len(witness["levels"]) == s and len(row["levels"]) == s, "incorrect level count")
    prod_a = prod_c = umax = 1
    for k, L in enumerate(witness["levels"], 1):
        require(L["k"] == k, "incorrect level order")
        theta, rho = (hex_integer(L[key]) for key in ("theta_hex", "rho_hex"))
        require(theta == Q // 2 ** k * prod_c and rho == prod_a, "level recursion failed")
        psi, determinant = reconstruct_adjugate(b, theta, rho, list(map(hex_integer, L["psi_seeds_hex"])))
        require(determinant == hex_integer(L["determinant_hex"]), "stored determinant differs")
        digest = hashlib.sha256()
        for x in psi:
            digest.update(f"{x:x}\n".encode("ascii"))
        require(digest.hexdigest() == L["sha256_psi"], "full adjugate hash mismatch")
        content = math.gcd(*psi)
        require(content == L["content"] == row["levels"][k-1]["cont"] == 1, "level is not primitive")
        u = psi[1:] + psi[:1]
        prod_a *= u[-1]
        prod_c *= 2 * u[0] + u[1]
        umax *= max(u)
    require(hashlib.sha256(str(umax).encode("ascii")).hexdigest() == row["sha256_max_u"], "max u hash mismatch")
    require(len(str(umax)) == row["digits_max_u"], "max u digit count differs")
    K = independent_K(b, s)
    require(K == Fraction(row["K"]), "K differs from the independently computed norms")
    require(Q >= 8 * K and Q - K >= 2 ** l, "explicit-set hypotheses fail")
    d = b ** s
    Delta = Fraction(2 ** b + 1, 2 ** b) ** ((d-1)//(b-1)) * Fraction(2,3) ** s
    require(Delta == Fraction(row["Delta_exact"]), "Delta differs")
    require(umax <= Q ** (d-1) * Delta * (1 + 10 * K / Q), "shape upper bound fails")
    exponent = l * (d-1)
    require(row["f_denominator_exponent"] == exponent, "incorrect f denominator")
    exact = Fraction(umax, 2 ** exponent)
    for value, key in ((exact, "f_bound"), (exact/2, "N_over_2^n")):
        require(re.fullmatch(r"[0-9]+\.[0-9]{6}", row[key]), "bound is not a six-place decimal string")
        gap = Fraction(row[key]) - value
        require(0 <= gap < Fraction(1, 10 ** 6), "bound is not the exact upward rounding")
    loss = Fraction(Q, 2 ** l) ** (d-1)
    require(0 <= Fraction(row["Q_loss_upper"]) - loss < Fraction(1, 10**9),
            "Q_loss_upper is not rounded upward to nine places")
    # (1+x)^m <= 1/(1-m*x), by binomial coefficients <= m^j.
    mx = (d-1) * Fraction(Q - 2 ** l, 2 ** l)
    require(0 <= mx < 1, "loss envelope is inapplicable")
    loss_envelope = 1 / (1 - mx)
    require(loss_envelope < Fraction(503,500), "Q loss is not certified below 1.006 < exp(.006)")
    return loss_envelope


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT / "results/certificates.json")
    args = parser.parse_args()
    rows = json.loads(args.manifest.read_text(encoding="utf-8"))
    expected = {(5,1),(7,1),(11,1),(13,1),(5,2),(7,2),(11,2),(13,2),
                (5,3),(7,3),(11,3),(13,3),(17,3),(13,4),(17,4)}
    require(len(rows) == 15 and {(r["b"],r["s"]) for r in rows} == expected, "the published 15 rows are incomplete")
    for row in rows:
        start = time.monotonic()
        envelope = verify_row(row, args.manifest.parent)
        print(f"PASS ({row['b']},{row['s']}): {row['s']} exact adjugates, f<={row['f_bound']}, "
              f"N/2^n<={row['N_over_2^n']}, loss<={envelope}<1.006 ({time.monotonic()-start:.2f}s)", flush=True)
    print("ALL PASS: 15 rows, 35 levels, independent determinants and K, full hashes, exact upper rounding.", flush=True)


if __name__ == "__main__":
    main()
