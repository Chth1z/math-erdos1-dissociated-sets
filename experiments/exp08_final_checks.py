"""Experiment 8: end-to-end brute-force confirmation for s=2, true dissociation threshold
for small cases, and certificates for larger (b,s)."""
import sys
from fractions import Fraction
from math import log2
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dissociated import is_Q_dissociated, is_dissociated_big as is_dissociated, lift, f_ratio  # noqa: E402
from tensor import K_bound, delta, tensor_certificate  # noqa: E402

sys.set_int_max_str_digits(0)


def brute_force_s2():
    print("== brute force: (b,s)=(5,2), l=1, |A|=25 ==")
    b, s = 5, 2
    K = K_bound(b, s)
    Q = 2 ** s
    found = 0
    while found < 3:
        Q += 2 ** (s + 1)
        if Q < 2 + K:
            continue
        cert = tensor_certificate(b, s, Q)
        if cert.primitive:
            u = cert.full_u()
            A = lift(u, 1)
            print(f"Q={Q}: primitive, |A|={len(A)}, dissociated={is_dissociated(A)}, f(A)={f_ratio(A):.4f}, "
                  f"min u={min(u)} max u={max(u)}")
            found += 1


def true_threshold():
    print("\n== true dissociation threshold vs guaranteed Q-K (b=5, s=1) ==")
    b = 5
    K = K_bound(b, 1)
    for Q in (16, 24, 40):
        cert = tensor_certificate(b, 1, Q)
        if not cert.primitive:
            continue
        u = cert.full_u()
        # largest Q' such that u is Q'-dissociated (brute force)
        Qp = 2
        while Qp <= 2 * Q and is_Q_dissociated(u, Qp):
            Qp += 1
        print(f"Q={Q}: u={u}, guaranteed (Q-K)={float(Q - K):.2f}, true threshold Q'={Qp - 1}")


def big_certificates():
    print("\n== certificates for larger (b,s) ==")
    for (b, s) in ((17, 3), (13, 4), (17, 4)):
        d = b ** s
        K = K_bound(b, s)
        Delta = delta(b, s)
        l = int(log2(200 * (d - 1) * float(K))) + 1
        Q0 = 2 ** (s + 1)
        t = int((2 ** l + K) // Q0) + 1
        tries = 0
        while True:
            tries += 1
            cert = tensor_certificate(b, s, Q0 * t)
            if cert.primitive:
                break
            t += 1
        umax = cert.u_max()
        f_exact = Fraction(umax, 2 ** (l * (d - 1)))
        print(f"(b,s)=({b},{s}) d={d} l={l} n={l*d}: Q={Q0*t} (tries={tries}) K={float(K):.1f} "
              f"Delta={float(Delta):.5f} f(n) <= {float(f_exact):.6f}  conts={[L.cont for L in cert.levels]}")


if __name__ == "__main__":
    brute_force_s2()
    true_threshold()
    big_certificates()
