"""Experiment 6: certificates.

(A) Theorem A check: rad(P_b) | q  =>  level-1 lattice primitive.  P_b = odd primes dividing
    b (2^b+1) prod_{0<i<j<b} N(1 + zeta^i + zeta^j).
(B) Coprimality certificates for levels k = 1..s-1 by counting: if gcd(a_k(t), c_k(t)) = 1 for
    more than 2*deg_t(a_k) integers t (with levels 1..k primitive), then a_k, c_k are coprime
    polynomials, hence gcd(a_k(t),c_k(t)) divides a fixed nonzero resultant for every t.
(C) Search small good Q and report a rational bound f(l d) <= max(u)/2^{l(d-1)}.
"""
import sys
from fractions import Fraction
from math import gcd
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tensor import K_bound, delta, tensor_certificate  # noqa: E402
from certificates import ceil_fraction, ceil_log2, upper_decimal  # noqa: E402


def rad_P(b):
    z = sp.symbols('z')
    Phi = sp.cyclotomic_poly(b, z)
    primes = set()
    for p in sp.factorint(b):
        primes.add(p)
    for p in sp.factorint(2 ** b + 1):
        primes.add(p)
    # norms N(1 + z^i + z^j) = resultant(Phi_b (or x^b-1 over all b-th roots), 1 + z^i + z^j)
    # use the full product over all b-th roots of unity != 1: resultant with (z^b-1)/(z-1)
    for i in range(1, b):
        for j in range(i + 1, b):
            r = sp.resultant(sp.Poly(sp.quo(z ** b - 1, z - 1), z), sp.Poly(1 + z ** i + z ** j, z))
            r = int(r)
            assert r != 0, "1+zeta^i+zeta^j vanishes: 3 | b?"
            for p in sp.factorint(abs(r)):
                primes.add(p)
    primes.discard(2)
    return sorted(primes)


def part_A():
    print("== (A) Theorem A: rad(P_b) | q forces level-1 primitivity ==")
    for b in (5, 7, 11):
        P = rad_P(b)
        M = 1
        for p in P:
            M *= p
        fails = 0
        tested = 0
        for m in range(1, 41):
            q = M * m
            cert = tensor_certificate(b, 1, 2 * q)
            tested += 1
            if not cert.primitive:
                fails += 1
        print(f"b={b}: P_b={P}, rad={M}; tested {tested} multiples q=rad*m: failures={fails}")


def level_degrees(b, s):
    """degrees in t of a_k(t), c_k(t) for Q = Q0*t: deg a_1 = b-1, deg theta_{k+1} = 1 + deg c_1 + ... + deg c_k,
    deg a_{k+1} = (b-1) * deg theta_{k+1} (homogeneous of degree b-1 in (theta, rho), deg rho < deg theta)."""
    degs = []
    deg_theta = 1
    for k in range(1, s + 1):
        d_ak = (b - 1) * deg_theta
        degs.append(d_ak)
        deg_theta = deg_theta + d_ak
    return degs


def part_B(b, s, Q0):
    print(f"\n== (B) coprimality certificate for (b,s)=({b},{s}), Q = {Q0}*t ==")
    degs = level_degrees(b, s)
    need = [2 * dk + 1 for dk in degs[:-1]]
    found = [0] * (s - 1)
    t = 1
    tried = 0
    while any(found[k] < need[k] for k in range(s - 1)) and tried < 20000:
        cert = tensor_certificate(b, s, Q0 * t)
        for k in range(s - 1):
            # levels 1..k+1 must be primitive for the polynomial values to be the actual ones
            if all(L.cont == 1 for L in cert.levels[:k + 1]):
                Lk = cert.levels[k]
                if gcd(Lk.a, Lk.c) == 1:
                    found[k] += 1
        t += 1
        tried += 1
    ok = all(found[k] >= need[k] for k in range(s - 1))
    print(f"   degrees of a_k: {degs[:-1]}, needed counts: {need}, found: {found}, "
          f"t values tried: {tried}  ->  coprimality certified: {ok}")
    return ok


def part_C(b, s, l, Q0, max_tries=5000):
    d = b ** s
    K = K_bound(b, s)
    Delta = delta(b, s)
    if l is None:
        # Initial target only; rounding/search may increase the final Q loss.
        l = ceil_log2(200 * (d - 1) * K)
    base = 2 ** l + K
    t0 = ceil_fraction(base / Q0)
    for t in range(t0, t0 + max_tries):
        Q = Q0 * t
        cert = tensor_certificate(b, s, Q)
        if cert.primitive:
            umax = cert.u_max()
            f_exact = Fraction(umax, 2 ** (l * (d - 1)))
            n = l * d
            print(f"(b,s)=({b},{s}) d={d} l={l} n={n}: Q={Q} (t={t}, tries={t-t0+1}) K={K} "
                  f"Delta~{float(Delta):.5f}  f(n) <= {upper_decimal(f_exact)}  "
                  f"[max u has {len(str(umax))} digits]")
            return cert, f_exact
    print(f"(b,s)=({b},{s}) l={l}: no primitive Q found in {max_tries} tries")
    return None, None


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    if "--no-AB" not in sys.argv:
        part_A()
        for (b, s) in ((5, 2), (7, 2), (5, 3), (11, 2), (11, 3), (13, 3)):
            Q0 = 2 ** (s + 1)
            part_B(b, s, Q0)
    print("\n== (C) explicit good Q and exact bounds ==")
    for (b, s) in ((5, 1), (7, 1), (11, 1), (5, 2), (7, 2), (11, 2), (7, 3), (11, 3), (13, 3)):
        part_C(b, s, None, 2 ** (s + 1))
