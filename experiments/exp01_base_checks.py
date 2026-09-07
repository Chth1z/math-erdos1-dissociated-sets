"""Experiment 1: sanity checks of the lattice construction and first primitivisation tests.

* verifies covol(Gamma) = h * Delta and the closed formulas from the exposition
* for the base pair (b = 3,5,7,9) and Q = 2^l, perturbs Q*Lambda by the identity block,
  reports gcd of minors, the ratio max(u)/Q^(d-1), and brute-force checks dissociation
  of the lifted set when it is small enough.
"""
import sys
from fractions import Fraction
from math import log2
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dissociated import f_ratio, is_dissociated, lift  # noqa: E402
from lattice import base_pair, choose_Q, dual_basis_l1, iterated_pair, primitivize  # noqa: E402


def check_pair(pair):
    h = pair.height()
    cg = pair.covol_gamma()
    d2 = pair.delta_sq()
    ok = (cg / h) ** 2 == d2
    print(f"{pair.label}: d={pair.d} den={pair.den} h={h} covol(Gamma)={cg} "
          f"Delta={float(cg / h):.6f}  [covol(Gamma)=h*Delta: {ok}]  sums zero: {pair.check_sum_zero()}")
    return ok


def main():
    print("== invariants of base and iterated pairs ==")
    for b in (3, 5, 7):
        check_pair(base_pair(b))
    for (b, s) in ((3, 2), (3, 3), (5, 2)):
        p = iterated_pair(b, s)
        check_pair(p)
        # closed formula: Delta = (1+2^-b)^{(b^s-1)/(b-1)} / (3/2)^s
        formula = Fraction(2 ** b + 1, 2 ** b) ** ((b ** s - 1) // (b - 1)) / Fraction(3, 2) ** s
        print(f"   closed formula Delta = {float(formula):.6f}")

    print("\n== primitivisation of the base pair ==")
    for b in (3, 5, 7, 9):
        pair = base_pair(b)
        norms, K = dual_basis_l1(pair)
        print(f"b={b}: dual basis l1 norms = {[str(x) for x in norms]}, K = {K} = {float(K):.3f}")
        for l in range(2, 13):
            Q = choose_Q(pair, l)
            P = primitivize(pair, Q)
            line = (f"   l={l:2d} Q={Q:6d} gcd={P.g} max(u)/Q^(d-1)={float(P.ratio):.6f} "
                    f"u={P.u if b <= 5 else '...'}")
            n = l * b
            if P.g == 1 and n <= 27:
                A = lift(P.u, l)
                ok = is_dissociated(A)
                line += f"  |A|={n} dissociated={ok} f(A)=max/2^(n-1)={f_ratio(A):.6f}"
            print(line)


if __name__ == "__main__":
    main()
