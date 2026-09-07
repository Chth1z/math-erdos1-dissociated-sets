"""Experiment 3: structure of the normal vector for iterated pairs, and sanity checks.

(a) For (b,s)=(5,2): is the normal vector u (identity-block perturbation) a tensor
    product u1 ⊗ u2 of two 5-vectors?  (Rank of u reshaped as 5x5.)
(b) Numerical admissibility check of Lambda_2 (b=3, d=9): no nonzero lattice vector in
    the open unit cube among small coefficient combinations.
(c) Base case: verify the polynomial criterion  gcd(minors)=1  <=>  for all primes p,
    deg gcd_p(phi, x^b-1) <= 1,  phi = q(2+x)(1-x)+1, on many (b,q).
"""
import itertools
import sys
from math import gcd
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lattice import base_pair, choose_Q, iterated_pair, primitivize  # noqa: E402


def part_a():
    print("== (a) tensor structure of u for (b,s)=(5,2) ==")
    pair = iterated_pair(5, 2)
    for l in (6, 8):
        Q = choose_Q(pair, l)
        P = primitivize(pair, Q)
        U = np.array(P.u, dtype=object).reshape(5, 5)
        M = sp.Matrix(U.tolist())
        print(f"l={l} Q={Q} gcd={P.g} rank(u as 5x5) = {M.rank()}")
        if M.rank() == 1:
            col = [int(x) for x in M[:, 0]]
            row = [int(x) for x in M[0, :]]
            g1 = 0
            for x in col:
                g1 = gcd(g1, x)
            g2 = 0
            for x in row:
                g2 = gcd(g2, x)
            print("   factor 1 (column, primitive):", [x // g1 for x in col])
            print("   factor 2 (row, primitive):   ", [x // g2 for x in row])
            print("   normalised ratios factor1:", [round(x / max(col), 5) for x in col])
            print("   normalised ratios factor2:", [round(x / max(row), 5) for x in row])


def part_b():
    print("\n== (b) admissibility sanity check of Lambda_2 (b=3, d=9) ==")
    pair = iterated_pair(3, 2)
    B = np.array(pair.basis, dtype=np.int64)          # 8 x 9 numerators, den 4
    den = pair.den
    # enumerate coefficient vectors m in [-M, M]^8 ; lattice vector = B^T m / den
    Mrange = 3
    rng = np.arange(-Mrange, Mrange + 1)
    grids = np.array(list(itertools.product(rng, repeat=B.shape[0])), dtype=np.int64)
    vecs = grids @ B                                  # numerators
    inside = np.all(np.abs(vecs) < den, axis=1)       # |x_i| < 1  <=>  |num| < den
    nonzero = np.any(grids != 0, axis=1)
    bad = inside & nonzero
    print(f"coefficient vectors tested: {len(grids)}; nonzero lattice vectors inside (-1,1)^9: {bad.sum()}")
    # real-t condition for v: t v in Lambda + (-1,1)^d  =>  |t| < 1 ; test t = 1, 1.25, 1.5
    v = np.array(pair.v, dtype=np.int64)
    for t_num in (4, 5, 6):                            # t = t_num/4
        # need: no lattice vector lam with |t v - lam|_inf < 1 (in numerators: |t_num v - 4 lam| < 4)
        target = t_num * v
        diff = np.abs(vecs - target)                   # vecs already numerators with den 4
        hit = np.all(diff < den, axis=1).sum()
        print(f"  t={t_num/4}: lattice points with t*v in lam + (-1,1)^d: {hit}  (must be 0 for t>=1)")


def poly_criterion(b, q):
    """True iff for all primes p, deg gcd_p(q(2+x)(1-x)+1, x^b-1) <= 1.
    Only primes dividing the resultant can be bad."""
    x = sp.symbols('x')
    phi = sp.Poly(q * (2 + x) * (1 - x) + 1, x)
    xb = sp.Poly(x**b - 1, x)
    res = int(sp.resultant(phi, xb))
    if res == 0:
        return False, []
    bad = []
    for p in sp.factorint(abs(res)):
        g = sp.gcd(sp.Poly(phi.as_expr(), x, modulus=p), sp.Poly(xb.as_expr(), x, modulus=p))
        if g.degree() >= 2:
            bad.append((p, g.degree()))
    return len(bad) == 0, bad


def part_c():
    print("\n== (c) base case: gcd of minors vs polynomial criterion ==")
    mismatches = 0
    for b in (5, 7, 9, 11, 13):
        pair = base_pair(b)
        for q in range(8, 40):
            Q = 2 * q
            P = primitivize(pair, Q)
            ok_poly, bad = poly_criterion(b, q)
            agree = (P.g == 1) == ok_poly
            if not agree:
                mismatches += 1
            if P.g != 1 or not ok_poly:
                print(f"b={b} q={q}: gcd={P.g}  poly-criterion ok={ok_poly} bad primes(deg)={bad}  agree={agree}")
    print(f"mismatches between gcd==1 and the polynomial criterion: {mismatches}")


if __name__ == "__main__":
    part_a()
    part_b()
    part_c()
