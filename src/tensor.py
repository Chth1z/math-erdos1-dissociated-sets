"""Tensor-structured primitivisation of the iterated lattices (polynomial model).

Ring R_b = Z[x]/(x^b - 1) ≅ Z^b (e_i <-> x^i).  Let g(x) = (2+x)(1-x) = 2 - x - x^2.

Level polynomials.  Given Q (even, with 2^s | Q) define recursively for k = 1..s
    theta_k = (Q / 2^k) * c_1 * ... * c_{k-1},     rho_k = a_1 * ... * a_{k-1},
    phi_k   = theta_k * g + rho_k,
    psi_k   = adj(phi_k~)  (adjugate of the multiplication-by-phi_k~ matrix, phi~(x)=phi(1/x)),
    u~_k    = x^{b-1} * psi_k / content(psi_k)        (primitive integer vector in Z^b),
    c_k     = 2 u~_{k,0} + u~_{k,1},   a_k = u~_{k,b-1}.
Then the perturbed lattice L ⊂ Z^{b^s} (basis  Q*lambda + monomial, see paper) has normal
vector  u = u~_1 ⊗ ... ⊗ u~_s,  and L is primitive iff content(psi_k) = 1 for all k.

When L is primitive the coordinates of u form a (Q - K)-dissociated set, where K bounds the
sum of the l1-norms of the dual basis:  K_1 exact,  K_s <= d'(1+4b/3) K_1 + (2b/3) K_{s-1}.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from functools import cache
from math import gcd
from typing import Sequence

import sympy as sp

from lattice import Pair, base_pair, bareiss_det, dual_basis_l1, minors_vector


# ----------------------------------------------------------------------------
# arithmetic in Z[x]/(x^b-1)
# ----------------------------------------------------------------------------

def cyc_mul(a: Sequence[int], c: Sequence[int]) -> list[int]:
    b = len(a)
    out = [0] * b
    for i, ai in enumerate(a):
        if ai:
            for j, cj in enumerate(c):
                if cj:
                    out[(i + j) % b] += ai * cj
    return out


def conj(a: Sequence[int]) -> list[int]:
    """a(x) -> a(1/x)."""
    b = len(a)
    return [a[(-i) % b] for i in range(b)]


def mult_matrix(a: Sequence[int]) -> list[list[int]]:
    """Matrix of multiplication by a on the basis 1, x, ..., x^{b-1} (columns = a*x^j)."""
    b = len(a)
    return [[a[(i - j) % b] for j in range(b)] for i in range(b)]


def content(v: Sequence[int]) -> int:
    g = 0
    for x in v:
        g = gcd(g, int(x))
    return g


def bareiss_solve_adj(M: Sequence[Sequence[int]], rhs: Sequence[int]) -> tuple[list[int], int]:
    """Fraction-free solution of  M x = det(M) * rhs  for an integer matrix M (so x = adj(M) rhs).

    Forward elimination is Bareiss (all intermediate entries are minors of [M|rhs]);
    back substitution uses only exact integer divisions.  Returns (x, det(M))."""
    n = len(M)
    A = [list(map(int, row)) + [int(rhs[i])] for i, row in enumerate(M)]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if swap is None:
                raise ZeroDivisionError("singular matrix")
            A[k], A[swap] = A[swap], A[k]
            sign = -sign
        piv = A[k][k]
        for i in range(k + 1, n):
            aik = A[i][k]
            row_i = A[i]
            row_k = A[k]
            for j in range(k + 1, n + 1):
                row_i[j] = (row_i[j] * piv - aik * row_k[j]) // prev
            row_i[k] = 0
        prev = piv
    D = sign * A[n - 1][n - 1]          # det(M)
    if D == 0:
        raise ZeroDivisionError("singular matrix")
    # U x = b'  with U upper triangular (rows of A), and x = adj(M) rhs / sign ... solve for y = |D| x
    # We solve U y = D_abs * b' with exact divisions, where D_abs = A[n-1][n-1] (= sign*D).
    Dabs = A[n - 1][n - 1]
    y = [0] * n
    for i in range(n - 1, -1, -1):
        num = Dabs * A[i][n]
        for j in range(i + 1, n):
            num -= A[i][j] * y[j]
        q, r = divmod(num, A[i][i])
        assert r == 0, "non-exact division in back substitution"
        y[i] = q
    # y = Dabs * M^{-1} rhs ;  adj(M) rhs = D * M^{-1} rhs = sign * y
    x = [sign * v for v in y]
    return x, D


def adjugate_element(a: Sequence[int]) -> tuple[list[int], int]:
    """psi with psi*a = N(a) (N = det of multiplication matrix); returns (psi, N).
    psi is the first column of adj(mult_matrix(a)); computed exactly with integers only."""
    b = len(a)
    M = mult_matrix(a)
    e0 = [1] + [0] * (b - 1)
    psi, N = bareiss_solve_adj(M, e0)
    if N == 0:
        raise ZeroDivisionError("a is a zero divisor")
    assert cyc_mul(psi, list(a)) == [N] + [0] * (b - 1)
    return psi, N


def g_poly(b: int) -> list[int]:
    g = [0] * b
    g[0] += 2
    g[1] -= 1
    g[2 % b] -= 1
    return g


@cache
def _shift_sum_charpoly(b: int) -> tuple[int, ...]:
    """Coefficients of det(z I - S - S^2), decreasing degree, for odd b.

    If r1,r2 solve x^2+x-z=0, T_j=r1^j+r2^j obeys T_j=-T_(j-1)+z*T_(j-2).
    The resultant with x^b-1 gives det(z I-S-S^2)=z^b+T_b(z)-1.
    """
    if b < 3 or b % 2 == 0:
        raise ValueError("require odd b >= 3")
    previous, current = [2], [-1]
    for _ in range(2, b + 1):
        following = [0] * max(len(current), len(previous) + 1)
        for i, x in enumerate(current):
            following[i] -= x
        for i, x in enumerate(previous):
            following[i + 1] += x
        previous, current = current, following
    coefficients = current + [0] * (b + 1 - len(current))
    coefficients[0] -= 1
    coefficients[b] += 1
    return tuple(reversed(coefficients))


def level_adjugate(b: int, theta: int, rho: int) -> tuple[list[int], int]:
    """Exact adj(m_conj(theta*g+rho)) e0 by Cayley-Hamilton, without elimination.

    Retains integers of the size of the output. adjugate_element remains the
    separate general-purpose Bareiss implementation for cross-checks.
    """
    if theta <= 0 or rho <= 0:
        raise ValueError("theta and rho must be positive")
    coefficients = _shift_sum_charpoly(b)
    a = 2 * theta + rho
    apowers, tpowers = [1], [1]
    for _ in range(b):
        apowers.append(apowers[-1] * a)
        tpowers.append(tpowers[-1] * theta)
    vector = [1] + [0] * (b - 1)
    psi = [apowers[b - 1]] + [0] * (b - 1)
    for k in range(1, b):
        vector = [vector[(i + 1) % b] + vector[(i + 2) % b] for i in range(b)]
        vector[0] += coefficients[k]
        weight = apowers[b - 1 - k] * tpowers[k]
        psi = [x + c * weight for x, c in zip(psi, vector)]
    determinant = sum(coefficients[k] * apowers[b - k] * tpowers[k] for k in range(b + 1))
    residual = [a * psi[i] - theta * psi[(i + 1) % b] - theta * psi[(i + 2) % b] for i in range(b)]
    if residual != [determinant] + [0] * (b - 1):
        raise ArithmeticError("the adjugate identity failed")
    return psi, determinant


# ----------------------------------------------------------------------------
# the recursion
# ----------------------------------------------------------------------------

@dataclass
class Level:
    k: int
    theta: int
    rho: int
    psi: list[int]
    N: int
    cont: int
    u: list[int]          # primitive vector x^{b-1} psi / cont, as coordinates u_0..u_{b-1}
    a: int
    c: int


@dataclass
class TensorCertificate:
    b: int
    s: int
    Q: int
    levels: list[Level] = field(default_factory=list)

    @property
    def primitive(self) -> bool:
        return all(L.cont == 1 for L in self.levels)

    def u_max(self) -> int:
        m = 1
        for L in self.levels:
            m *= max(L.u)
        return m

    def u_min(self) -> int:
        m = 1
        for L in self.levels:
            m *= min(L.u)
        return m

    def full_u(self) -> list[int]:
        """Kronecker product u~_1 ⊗ ... ⊗ u~_s (length b^s).  Only for small cases."""
        vec = [1]
        for L in self.levels:
            vec = [x * y for x in vec for y in L.u]
        return vec


def tensor_certificate(b: int, s: int, Q: int, *, method: str = "charpoly") -> TensorCertificate:
    if b < 3 or b % 2 == 0 or s < 1 or Q <= 0 or Q % (2 ** s):
        raise ValueError("require odd b >= 3, s >= 1, and positive Q divisible by 2^s")
    if method not in ("charpoly", "bareiss"):
        raise ValueError("method must be 'charpoly' or 'bareiss'")
    g = g_poly(b)
    cert = TensorCertificate(b, s, Q)
    prod_c, prod_a = 1, 1
    for k in range(1, s + 1):
        theta = (Q // 2 ** k) * prod_c
        rho = prod_a
        phi = [theta * gi for gi in g]
        phi[0] += rho
        psi, N = level_adjugate(b, theta, rho) if method == "charpoly" else adjugate_element(conj(phi))
        ct = content(psi)
        # u = x^{b-1} psi / ct  :  coordinate i of x^{b-1} psi is psi_{i+1 mod b}
        u = [psi[(i + 1) % b] // ct for i in range(b)]
        if all(x <= 0 for x in u):
            u = [-x for x in u]
        if not all(x > 0 for x in u):
            raise ArithmeticError("the level normal does not have positive coordinates")
        a = u[b - 1]
        c = 2 * u[0] + u[1]
        cert.levels.append(Level(k, theta, rho, psi, N, ct, u, a, c))
        prod_c *= c
        prod_a *= a
    return cert


# ----------------------------------------------------------------------------
# the perturbed lattice itself (for cross-checks in small dimension)
# ----------------------------------------------------------------------------

def monomial_index(exps: Sequence[int], b: int) -> int:
    """x_1^{e_1} ... x_s^{e_s}  ->  index with x_1 slowest (matches Kronecker order of full_u)."""
    idx = 0
    for e in exps:
        idx = idx * b + e
    return idx


def poly_to_vec(coeffs: dict[tuple[int, ...], Fraction], b: int, s: int) -> list[Fraction]:
    v = [Fraction(0)] * (b ** s)
    for exps, cval in coeffs.items():
        v[monomial_index(exps, b)] += cval
    return v


def mul_multivar(p: dict, q: dict, b: int) -> dict:
    out: dict[tuple[int, ...], Fraction] = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            e = tuple((x + y) % b for x, y in zip(e1, e2))
            out[e] = out.get(e, Fraction(0)) + c1 * c2
    return {e: c for e, c in out.items() if c != 0}


def perturbed_basis(b: int, s: int, Q: int) -> tuple[list[list[int]], list[list[int]]]:
    """Return (unperturbed basis W of Q*Lambda_s, perturbed basis W+E) as integer column vectors
    in Z^{b^s}, following the monomial matching of the paper."""
    import itertools

    def mono(var: int, e: int) -> dict:
        exps = [0] * s
        exps[var] = e % b
        return {tuple(exps): Fraction(1)}

    def const(cval) -> dict:
        return {tuple([0] * s): Fraction(cval)}

    def add(p: dict, q: dict) -> dict:
        out = dict(p)
        for e, c in q.items():
            out[e] = out.get(e, Fraction(0)) + c
        return {e: c for e, c in out.items() if c != 0}

    W, WE = [], []
    for k in range(1, s + 1):
        # factor  Q 2^{-k} (2+x_1)...(2+x_k)(1-x_k)
        f = const(Fraction(Q, 2 ** k))
        for j in range(k):
            f = mul_multivar(f, add(const(2), mono(j, 1)), b)
        f = mul_multivar(f, add(const(1), {tuple([0] * (k - 1) + [1] + [0] * (s - k)): Fraction(-1)}), b)
        # monomials x_k^i m(x_{>k}),  i <= b-2
        for i in range(b - 1):
            for rest in itertools.product(range(b), repeat=s - k):
                exps = [0] * s
                exps[k - 1] = i
                for t, e in enumerate(rest):
                    exps[k + t] = e
                m = {tuple(exps): Fraction(1)}
                vec = poly_to_vec(mul_multivar(f, m, b), b, s)
                assert all(x.denominator == 1 for x in vec)
                col = [int(x) for x in vec]
                W.append(col)
                pert_exps = [b - 1] * (k - 1) + exps[k - 1:]
                col2 = list(col)
                col2[monomial_index(pert_exps, b)] += 1
                WE.append(col2)
    return W, WE


# ----------------------------------------------------------------------------
# constants
# ----------------------------------------------------------------------------

def _base_dual_vectors(b: int) -> list[list[Fraction]]:
    """Dual basis vectors b_i^* in V_b of the base basis b_i = (1/2)(2+x)(1-x)x^i, i <= b-2."""
    from lattice import solve_rational
    pair = base_pair(b)
    B = [[Fraction(x, pair.den) for x in col] for col in pair.basis]
    n = len(B)
    G = [[sum(a * c for a, c in zip(B[i], B[j])) for j in range(n)] for i in range(n)]
    duals = []
    for j in range(n):
        e = [Fraction(int(i == j)) for i in range(n)]
        coef = solve_rational(G, e)
        duals.append([sum(coef[i] * B[i][k] for i in range(n)) for k in range(b)])
    return duals


def K_bound(b: int, s: int) -> Fraction:
    """Upper bound K_s for the sum of the l1 norms of a dual basis of the tensor basis of
    Lambda_s.  Recursion (exact for the explicit dual vectors constructed in the paper):
        first family  (b_i ⊗ e_m)^* = b_i^* ⊗ e_m - c_i 1_b ⊗ (e_m - 1_{d'}/d'),  c_i = <b_i^*, v_1>/h_1,
        second family (v_1 ⊗ b'_j)^* = (1_b/h_1) ⊗ (b'_j)^*.
    """
    duals = _base_dual_vectors(b)
    h1 = Fraction(3, 2)
    v1 = [Fraction(0)] * b
    v1[0], v1[1] = Fraction(1), Fraction(1, 2)
    K1 = sum(sum(abs(x) for x in D) for D in duals)
    K = K1
    for level in range(2, s + 1):
        dprime = b ** (level - 1)
        first = Fraction(0)
        for D in duals:
            c = sum(x * y for x, y in zip(D, v1)) / h1
            shift = c * (1 - Fraction(1, dprime))
            block_m = sum(abs(x - shift) for x in D)
            others = abs(c) * b * Fraction(dprime - 1, dprime)
            first += block_m + others
        K = dprime * first + (b / h1) * K
    return K


def delta(b: int, s: int) -> Fraction:
    return Fraction(2 ** b + 1, 2 ** b) ** ((b ** s - 1) // (b - 1)) / Fraction(3, 2) ** s
