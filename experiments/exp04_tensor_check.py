"""Experiment 4: cross-check the tensor recursion against direct minors.

For (b,s) in {(5,1),(5,2),(7,2),(3,2)} and several Q:
  * build the perturbed basis W+E of Q*Lambda_s with the monomial matching
  * compute its minors vector directly (gcd, direction)
  * compare with the tensor certificate: gcd == prod content(psi_k), direction == ⊗ u~_k
Also: the unperturbed W must have all minors equal (normal ∝ 1), i.e. W spans Q*Lambda_s ⊂ V.
"""
import sys
from math import gcd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dissociated import is_dissociated, lift, f_ratio  # noqa: E402
from lattice import minors_vector  # noqa: E402
from tensor import perturbed_basis, tensor_certificate, K_bound, delta  # noqa: E402


def main():
    for (b, s) in ((5, 1), (5, 2), (3, 2), (7, 2)):
        for Q in (2 ** s * 3, 2 ** s * 5, 2 ** s * 8, 2 ** s * 11):
            W, WE = perturbed_basis(b, s, Q)
            d = b ** s
            assert len(W) == d - 1 and all(sum(c) == 0 for c in W)
            w0 = minors_vector(W)
            assert len(set(abs(x) for x in w0)) == 1, "unperturbed minors not all equal"
            w = minors_vector(WE)
            g = 0
            for x in w:
                g = gcd(g, x)
            wp = [x // g for x in w]
            if all(x <= 0 for x in wp):
                wp = [-x for x in wp]
            cert = tensor_certificate(b, s, Q)
            u = cert.full_u()
            prod_cont = 1
            for L in cert.levels:
                prod_cont *= L.cont
            same_dir = (wp == u)
            print(f"(b,s)=({b},{s}) Q={Q:3d}: direct gcd={g:6d}  prod content={prod_cont:6d}  "
                  f"direction match={same_dir}  conts={[L.cont for L in cert.levels]}  "
                  f"max(u)/Q^(d-1)={cert.u_max()/Q**(d-1):.5f}  Delta={float(delta(b,s)):.5f}")
    print("\nK bounds:", {(b, s): float(K_bound(b, s)) for (b, s) in ((5, 1), (5, 2), (7, 2), (11, 3), (13, 3), (11, 4))})

    print("\n== small explicit dissociated sets from the tensor construction ==")
    for (b, s, l) in ((5, 1, 3), (5, 1, 4), (5, 1, 5), (7, 1, 3), (5, 2, 1)):
        K = K_bound(b, s)
        Q = 2 ** l + int(K) + 1
        Q += (-Q) % (2 ** s)
        cert = tensor_certificate(b, s, Q)
        if not cert.primitive:
            print(f"(b,s,l)=({b},{s},{l}) Q={Q}: not primitive, conts={[L.cont for L in cert.levels]}")
            continue
        u = cert.full_u()
        A = lift(u, l)
        ok = is_dissociated(A) if len(A) <= 27 else None
        print(f"(b,s,l)=({b},{s},{l}) Q={Q} |A|={len(A)} dissociated={ok} f(A)={f_ratio(A):.5f} "
              f"max(u)/Q^(d-1)={cert.u_max()/Q**(b**s-1):.5f}")


if __name__ == "__main__":
    main()
