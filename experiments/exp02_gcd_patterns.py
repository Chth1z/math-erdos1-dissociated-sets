"""Experiment 2: when is the identity-block perturbation primitive (gcd = 1)?

* base pairs b = 3..21 (odd), several Q
* alternative perturbation rows for the failing cases
* iterated pairs (b,s) = (5,2), (7,2), (3,2), (3,3)
Also records K (dual basis constant) to see how it grows.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lattice import base_pair, choose_Q, dual_basis_l1, iterated_pair, primitivize  # noqa: E402


def gcds_for(pair, ls, rows=None):
    out = []
    for l in ls:
        Q = choose_Q(pair, l)
        P = primitivize(pair, Q, rows)
        out.append((l, Q, P.g, float(P.ratio)))
    return out


def main():
    print("== base pairs, identity block on rows 0..d-2 ==")
    for b in range(3, 23, 2):
        pair = base_pair(b)
        _, K = dual_basis_l1(pair)
        res = gcds_for(pair, (4, 6, 8, 10))
        gs = [g for (_, _, g, _) in res]
        print(f"b={b:2d} K={float(K):7.3f} Delta={float(pair.covol_gamma()/pair.height()):.5f} gcds={gs}"
              f"  ratio@l=10={res[-1][3]:.5f}")

    print("\n== failing b: try other perturbation-row patterns ==")
    d = 9
    patterns = {
        "rows 1..d-1 (skip row 0)": list(range(1, d)),
        "rows shifted by 2": [(j + 2) % d for j in range(d - 1)],
        "reverse rows d-1..1": [d - 1 - j for j in range(d - 1)],
        "row j -> (2j) mod d": [(2 * j) % d for j in range(d - 1)],
        "row j -> (4j) mod d": [(4 * j) % d for j in range(d - 1)],
    }
    for b in (3, 9, 15, 21):
        pair = base_pair(b)
        d = b
        pats = {
            "rows 0..d-2": list(range(d - 1)),
            "rows 1..d-1": list(range(1, d)),
            "rows shifted by 2": [(j + 2) % d for j in range(d - 1)],
            "reverse": [d - 1 - j for j in range(d - 1)],
            "j->2j mod d": [(2 * j) % d for j in range(d - 1)],
            "j->4j mod d": [(4 * j) % d for j in range(d - 1)],
            "j->(j+1)": [(j + 1) % d for j in range(d - 1)],
        }
        for name, rows in pats.items():
            if len(set(rows)) != d - 1:
                continue
            res = gcds_for(pair, (4, 6, 8), rows)
            print(f"b={b:2d} {name:20s} gcds={[g for (_, _, g, _) in res]}")

    print("\n== iterated pairs ==")
    for (b, s) in ((3, 2), (5, 2), (7, 2), (3, 3)):
        pair = iterated_pair(b, s)
        _, K = dual_basis_l1(pair)
        res = gcds_for(pair, (6, 8, 10))
        print(f"(b,s)=({b},{s}) d={pair.d} den={pair.den} K={float(K):8.3f} "
              f"Delta={float(pair.covol_gamma()/pair.height()):.5f} "
              f"gcds={[g for (_, _, g, _) in res]} ratios={[round(r, 5) for (_, _, _, r) in res]}")


if __name__ == "__main__":
    main()
