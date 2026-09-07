"""Experiment 5: the level polynomials a_k(t), c_k(t) as polynomials in the free parameter.

Level 1 with Q = 2q:  phi_1 = q g + 1,  psi_1 = adj(phi_1~) has entries in Z[q];
    a_1(q) = psi_{1,0},  c_1(q) = 2 psi_{1,1} + psi_{1,2}.
We check: gcd(a_1, c_1) in Q[q], factorisation of a_1, the identity a_1 = N - (q/b) N',
and the resultant Res_q(a_1, c_1) (its prime factors are the only primes that can ever
divide gcd(a_1(q), c_1(q))).
Level 2 (b=5,7) with Q = 8t: a_2(t), c_2(t) and their gcd / resultant.
"""
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tensor import g_poly  # noqa: E402

x, q, t = sp.symbols('x q t')


def circulant_conj_matrix(phi_coeffs, b):
    """multiplication matrix of phi~(x) = phi(1/x) in Z[x]/(x^b-1); phi given as list of sympy exprs."""
    conj = [phi_coeffs[(-i) % b] for i in range(b)]
    return sp.Matrix(b, b, lambda i, j: conj[(i - j) % b])


def level_polys(b, theta, rho):
    """Return (psi list of sympy exprs, N) for phi = theta*g + rho."""
    g = g_poly(b)
    phi = [sp.expand(theta * gi) for gi in g]
    phi[0] = sp.expand(phi[0] + rho)
    M = circulant_conj_matrix(phi, b)
    N = sp.factor(M.det(method="berkowitz"))
    adj = M.adjugate()
    psi = [sp.expand(adj[i, 0]) for i in range(b)]
    return psi, sp.expand(N)


def analyse_level1(b):
    psi, N = level_polys(b, q, 1)
    a1 = psi[0]
    c1 = sp.expand(2 * psi[1] + psi[2])
    gcd_poly = sp.gcd(sp.Poly(a1, q), sp.Poly(c1, q))
    res = sp.resultant(sp.Poly(a1, q), sp.Poly(c1, q))
    ident = sp.simplify(a1 - (N - q * sp.diff(N, q) / b))
    print(f"b={b}: deg a1={sp.degree(a1, q)}, deg c1={sp.degree(c1, q)}, gcd={gcd_poly.as_expr()}, "
          f"a1 = N - (q/b)N' : {ident == 0}")
    print(f"     a1 factors: {sp.factor(a1)}")
    print(f"     c1 - 3 a1 factors: {sp.factor(sp.expand(c1 - 3 * a1))}")
    print(f"     Res(a1,c1) = {res} = {sp.factorint(int(res)) if abs(int(res)) < 10**40 else 'big'}")
    return a1, c1, psi


def analyse_level2(b):
    psi1, _ = level_polys(b, 2 * t, 1)          # Q = 4t*? we use Q = 8t: theta_1 = Q/2 = 4t
    # careful: theta_1 = Q/2.  Take Q = 8t -> theta_1 = 4t.
    psi1, _ = level_polys(b, 4 * t, 1)
    a1 = psi1[0]
    c1 = sp.expand(2 * psi1[1] + psi1[2])
    # content of psi1 as polynomial vector must be 1 generically; assume primitive level 1
    theta2 = sp.expand(2 * t * c1)              # Q/4 * c1 = 2t c1
    rho2 = a1
    psi2, _ = level_polys(b, theta2, rho2)
    a2 = psi2[0]
    c2 = sp.expand(2 * psi2[1] + psi2[2])
    P_a2, P_c2 = sp.Poly(a2, t), sp.Poly(c2, t)
    gcd_poly = sp.gcd(P_a2, P_c2)
    print(f"b={b} level 2: deg a2={P_a2.degree()}, deg c2={P_c2.degree()}, gcd = {gcd_poly.as_expr()}")
    res = sp.resultant(P_a2, P_c2)
    r = int(res)
    print(f"     Res(a2,c2) has {len(str(abs(r)))} digits; nonzero: {r != 0}")
    # small prime factors of the resultant (trial division up to 10^5)
    small = []
    rr = abs(r)
    for p in sp.primerange(2, 100000):
        if rr % p == 0:
            e = 0
            while rr % p == 0:
                rr //= p
                e += 1
            small.append((p, e))
    print(f"     small prime factors (<1e5) of Res: {small}")
    return a2, c2


if __name__ == "__main__":
    for b in (5, 7, 11, 13):
        analyse_level1(b)
    for b in (5, 7):
        analyse_level2(b)
