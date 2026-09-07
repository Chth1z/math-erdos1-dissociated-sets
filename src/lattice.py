"""Explicit version of the lattice construction behind the disproof of Erdős #1.

Notation follows T. Bloom's exposition (erdosproblems.com/1, Sept 2026) of the
GPT-6 Astra construction:

* V_d = {x in R^d : sum x_i = 0}.
* An *admissible pair* (Lambda, v): Lambda a full-rank lattice in V_d with
  Lambda ∩ (-1,1)^d = {0}, and v ∉ V_d such that t v ∈ Lambda + (-1,1)^d implies |t|<1.
  Height h = sum_i v_i.  Gamma = Lambda + Z v has covol(Gamma) = h * Delta where
  Delta = covol(Lambda)/sqrt(d).
* Base pair in R^b (b odd): Lambda_1 = M(Z^b ∩ V_b), v_1 = M e_0, M = I + P/2,
  P the cyclic shift (P e_i = e_{i+1}).  h_1 = 3/2, covol(Gamma_1) = 1 + 2^-b.
* Iteration: from (Lambda, v) in R^d build a pair in R^{bd}:
      Lambda' = {(lambda_0 + beta_0 v, ..., lambda_{b-1} + beta_{b-1} v)},
      v'      = v_1 ⊗ v,
  heights multiply, covol(Gamma') = covol(Gamma_1) covol(Gamma)^b.

All vectors are stored as integer numerator lists with a common power-of-two
denominator `den`, so every computation below is exact.

The new ingredient here (not in the exposition) is `primitivize`, which turns
Q*Lambda into an explicit *primitive* integer lattice by adding the identity block
to a basis matrix, and returns the integer normal vector u whose coordinates form a
(Q-K)-dissociated set, K an explicit constant depending only on Lambda.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd, isqrt
from typing import Sequence

Vec = list[int]


# ----------------------------------------------------------------------------
# exact linear algebra on Python integers
# ----------------------------------------------------------------------------

def bareiss_det(M: Sequence[Sequence[int]]) -> int:
    """Fraction-free Gaussian elimination; exact determinant of an integer matrix."""
    n = len(M)
    if n == 0:
        return 1
    A = [list(map(int, row)) for row in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if swap is None:
                return 0
            A[k], A[swap] = A[swap], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
            A[i][k] = 0
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def gram_det(cols: Sequence[Vec]) -> int:
    """det(B^T B) for integer column vectors (exact)."""
    G = [[sum(a * b for a, b in zip(c1, c2)) for c2 in cols] for c1 in cols]
    return bareiss_det(G)


def minors_vector(cols: Sequence[Vec]) -> Vec:
    """For a d x (d-1) integer matrix (given by columns) return the vector w with
    w_i = (-1)^i det(matrix with row i deleted).  w is orthogonal to every column
    and gcd(w) = index of the column span in its saturation."""
    d = len(cols[0])
    assert len(cols) == d - 1
    rows = [[c[i] for c in cols] for i in range(d)]
    w = []
    for i in range(d):
        sub = rows[:i] + rows[i + 1:]
        w.append((-1) ** i * bareiss_det(sub))
    return w


def solve_rational(A: Sequence[Sequence[Fraction]], b: Sequence[Fraction]) -> list[Fraction]:
    """Solve A x = b exactly (A square, nonsingular) by Gaussian elimination."""
    n = len(A)
    M = [list(row) + [b[i]] for i, row in enumerate(A)]
    for k in range(n):
        piv = next(i for i in range(k, n) if M[i][k] != 0)
        M[k], M[piv] = M[piv], M[k]
        inv = 1 / M[k][k]
        M[k] = [x * inv for x in M[k]]
        for i in range(n):
            if i != k and M[i][k] != 0:
                f = M[i][k]
                M[i] = [x - f * y for x, y in zip(M[i], M[k])]
    return [M[i][n] for i in range(n)]


# ----------------------------------------------------------------------------
# admissible pairs
# ----------------------------------------------------------------------------

@dataclass
class Pair:
    d: int              # ambient dimension
    den: int            # common denominator (power of two)
    basis: list[Vec]    # numerators of a basis of Lambda (d-1 vectors in V_d)
    v: Vec              # numerators of v
    label: str = ""

    # -- basic invariants -------------------------------------------------
    def height(self) -> Fraction:
        return Fraction(sum(self.v), self.den)

    def covol_lambda_sq(self) -> Fraction:
        return Fraction(gram_det(self.basis), self.den ** (2 * (self.d - 1)))

    def delta_sq(self) -> Fraction:
        """Delta^2 = covol(Lambda)^2 / d."""
        return self.covol_lambda_sq() / self.d

    def covol_gamma(self) -> Fraction:
        cols = self.basis + [self.v]
        rows = [[c[i] for c in cols] for i in range(self.d)]
        return Fraction(abs(bareiss_det(rows)), self.den ** self.d)

    def check_sum_zero(self) -> bool:
        return all(sum(c) == 0 for c in self.basis)


def base_pair(b: int) -> Pair:
    """Bloom's base pair in R^b, b odd: basis 2*M(e_i - e_{i+1}) = 2e_i - e_{i+1} - e_{i+2},
    2*v_1 = 2e_0 + e_1, denominator 2."""
    assert b % 2 == 1 and b >= 3
    basis = []
    for i in range(b - 1):
        c = [0] * b
        c[i] += 2
        c[(i + 1) % b] -= 1
        c[(i + 2) % b] -= 1
        basis.append(c)
    v = [0] * b
    v[0] = 2
    v[1] = 1
    return Pair(b, 2, basis, v, label=f"base(b={b})")


def iterate(pair: Pair, base: Pair) -> Pair:
    """One step of the tensor iteration: (Lambda, v) in R^d, base (Lambda_1, v_1) in R^b
    -> pair in R^{bd}.  Denominators multiply."""
    d, b = pair.d, base.d
    D, D1 = pair.den, base.den
    newden = D * D1
    basis: list[Vec] = []
    # block copies of Lambda (scaled to the new denominator)
    for j in range(b):
        for lam in pair.basis:
            c = [0] * (b * d)
            for i, x in enumerate(lam):
                c[j * d + i] = x * D1
            basis.append(c)
    # coupling vectors beta ⊗ v for beta in basis of Lambda_1
    for beta in base.basis:
        c = [0] * (b * d)
        for j in range(b):
            for i, x in enumerate(pair.v):
                c[j * d + i] = beta[j] * x
        basis.append(c)
    v = [0] * (b * d)
    for j in range(b):
        for i, x in enumerate(pair.v):
            v[j * d + i] = base.v[j] * x
    return Pair(b * d, newden, basis, v, label=f"{pair.label}⊗{base.label}")


def iterated_pair(b: int, s: int) -> Pair:
    base = base_pair(b)
    pair = base
    for _ in range(s - 1):
        pair = iterate(pair, base)
    return pair


# ----------------------------------------------------------------------------
# dual basis and the perturbation constant K
# ----------------------------------------------------------------------------

def dual_basis_l1(pair: Pair) -> tuple[list[Fraction], Fraction]:
    """ℓ1 norms of the dual basis vectors b_j^* ∈ V_d (⟨b_i, b_j^*⟩ = δ_ij) and their sum K.

    For any λ = Σ m_j b_j ∈ Λ we have |m_j| ≤ |b_j^*|_1 |λ|_∞, hence the lattice
    span(Q b_j + e_j) with |e_j|_∞ ≤ 1 avoids (-(Q-K), Q-K)^d whenever Q > K."""
    B = [[Fraction(x, pair.den) for x in col] for col in pair.basis]   # columns
    n = len(B)
    G = [[sum(a * c for a, c in zip(B[i], B[j])) for j in range(n)] for i in range(n)]
    # dual basis: b_j^* = Σ_i (G^{-1})_{ij} b_i ; compute G^{-1} column by column
    norms = []
    for j in range(n):
        e = [Fraction(int(i == j)) for i in range(n)]
        coef = solve_rational(G, e)
        bstar = [sum(coef[i] * B[i][k] for i in range(n)) for k in range(pair.d)]
        norms.append(sum(abs(x) for x in bstar))
    return norms, sum(norms)


# ----------------------------------------------------------------------------
# primitivisation: explicit integer normal vector
# ----------------------------------------------------------------------------

@dataclass
class Primitive:
    Q: int
    u: Vec              # integer normal vector (all coordinates positive)
    g: int              # gcd of minors (1 means the perturbed lattice is primitive)
    K: Fraction         # perturbation constant; u is (Q-K)-dissociated when g == 1
    ratio: Fraction     # max(u) / Q^(d-1)


def primitivize(pair: Pair, Q: int, perturb_rows: Sequence[int] | None = None) -> Primitive:
    """Perturb the basis Q*b_j by adding e_{r_j} (r_j = perturb_rows[j], default j) and
    return the minors vector of the resulting d x (d-1) matrix.

    Requires den | Q.  If the returned g equals 1, the coordinates of u form a
    (Q-K)-dissociated set of positive integers, K = dual_basis_l1(pair)[1]."""
    assert Q % pair.den == 0
    d = pair.d
    scale = Q // pair.den
    cols = [[x * scale for x in c] for c in pair.basis]
    rows_to_perturb = list(perturb_rows) if perturb_rows is not None else list(range(d - 1))
    for j, r in enumerate(rows_to_perturb):
        cols[j][r] += 1
    w = minors_vector(cols)
    g = 0
    for x in w:
        g = gcd(g, x)
    if all(x <= 0 for x in w):
        w = [-x for x in w]
    assert all(x > 0 for x in w), "minors vector not of constant sign; Q too small?"
    _, K = dual_basis_l1(pair)
    ratio = Fraction(max(w), Q ** (d - 1))
    return Primitive(Q, w, g, K, ratio)


def choose_Q(pair: Pair, l: int) -> int:
    """Smallest multiple of den that is >= 2^l + K, so that u is 2^l-dissociated and the
    lift with l doublings is a genuine dissociated set of size l*d."""
    _, K = dual_basis_l1(pair)
    need = 2 ** l + K
    q = -(-need // pair.den) * pair.den  # ceil to a multiple of den
    return int(q)
