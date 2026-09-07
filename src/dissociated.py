"""Utilities for dissociated (distinct-subset-sum) sets.

A finite set A of positive integers is *dissociated* if all 2^|A| subset sums
are distinct.  Erdős' problem #1 asks whether max A >> 2^|A| must hold.

Functions
---------
is_dissociated(A)            exact check (meet in the middle on signed sums)
is_Q_dissociated(u, Q)       exact check that no nonzero c in (-Q,Q)^d has c.u = 0
                             (brute force, small d*Q only)
lift(U, l)                   {2^j * u : 0 <= j < l, u in U}; dissociated whenever
                             U is 2^l-dissociated
f_ratio(A)                   max(A) / 2^(|A|-1)
"""
from __future__ import annotations

import itertools
from math import gcd
from typing import Iterable, Sequence

import numpy as np


def signed_sums(a: Sequence[int]) -> np.ndarray:
    """All 3^n signed sums sum(eps_i a_i), eps_i in {-1,0,1}, as an int64 array."""
    s = np.zeros(1, dtype=np.int64)
    for x in a:
        s = np.concatenate((s - x, s, s + x))
    return s


def is_dissociated(A: Iterable[int]) -> bool:
    """Exact test.  A is dissociated iff no nontrivial signed sum vanishes.

    Splits A into two halves; the signed sums of each half have 3^(n/2) elements,
    so n <= 30 or so is comfortable.  Values must fit in int64 (max |sum| < 2^63).
    """
    A = sorted(int(x) for x in A)
    if len(A) != len(set(A)) or any(x <= 0 for x in A):
        return False
    if len(A) <= 1:
        return True
    if sum(A) >= 2**62:
        raise ValueError("elements too large for int64 arithmetic; use is_dissociated_big")
    h = len(A) // 2
    d1 = signed_sums(A[:h])
    d2 = signed_sums(A[h:])
    # each half must itself be dissociated: zero attained only by the trivial signing
    if np.count_nonzero(d1 == 0) != 1 or np.count_nonzero(d2 == 0) != 1:
        return False
    common = np.intersect1d(d1, -d2)
    return common.size == 1 and common[0] == 0


def is_dissociated_big(A: Iterable[int]) -> bool:
    """Exact meet-in-the-middle test with Python integers (arbitrary size), n <= ~28.

    A is dissociated iff the only signed sum  sum eps_i a_i  (eps_i in {-1,0,1}) that
    vanishes is the trivial one.  Split A = A1 ∪ A2; enumerate signed sums of each half
    as multisets; require 0 to occur once in each and D1 ∩ (-D2) = {0}."""
    A = sorted(int(x) for x in A)
    if len(A) != len(set(A)) or any(x <= 0 for x in A):
        return False
    if len(A) <= 1:
        return True
    h = len(A) // 2

    def signed(part):
        sums = [0]
        for x in part:
            sums = [s - x for s in sums] + sums + [s + x for s in sums]
        return sums

    d1 = signed(A[:h])
    d2 = signed(A[h:])
    if d1.count(0) != 1 or d2.count(0) != 1:
        return False
    s1 = set(d1)
    for v in d2:
        if -v in s1 and v != 0:
            return False
    return True


def is_Q_dissociated(u: Sequence[int], Q: int) -> bool:
    """Brute force: no nonzero integer vector c with |c_i| < Q and sum c_i u_i = 0.

    Enumerates (2Q-1)^(d-1) vectors for the first d-1 coordinates and solves for the
    last; only for tiny instances.
    """
    u = [int(x) for x in u]
    d = len(u)
    rng = range(-(Q - 1), Q)
    for c in itertools.product(rng, repeat=d - 1):
        s = sum(ci * ui for ci, ui in zip(c, u[:-1]))
        # need c_d * u_d = -s with |c_d| < Q
        if s % u[-1] == 0:
            cd = -s // u[-1]
            if abs(cd) < Q and (cd != 0 or any(c)):
                return False
    return True


def lift(U: Sequence[int], l: int) -> list[int]:
    """A = {2^j u : 0<=j<l, u in U}.  If U is 2^l-dissociated then A is dissociated
    with |A| = l*|U| and max A = 2^(l-1) max U."""
    return sorted((1 << j) * int(u) for u in U for j in range(l))


def f_ratio(A: Sequence[int]) -> float:
    """max(A) / 2^(|A|-1); Erdős' problem asks whether this is bounded below."""
    A = list(A)
    return max(A) / 2 ** (len(A) - 1)


def vec_gcd(v: Iterable[int]) -> int:
    g = 0
    for x in v:
        g = gcd(g, int(x))
    return g
