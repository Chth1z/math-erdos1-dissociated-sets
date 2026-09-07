"""Experiment 7: how loose is the recursive bound K_s for the tensor basis?
Compute the exact sum of l1-norms of the dual basis for (b,s) = (5,2), (7,2), (5,3)
and compare with K_bound."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lattice import Pair, dual_basis_l1  # noqa: E402
from tensor import K_bound, perturbed_basis  # noqa: E402


def exact_K(b, s):
    Q = 2 ** s  # W = Q * basis with den 2^s  ->  Pair with den = Q gives the true basis
    W, _ = perturbed_basis(b, s, Q)
    d = b ** s
    v = [0] * d  # unused
    pair = Pair(d, Q, W, v)
    norms, K = dual_basis_l1(pair)
    return K, max(norms)


if __name__ == "__main__":
    for (b, s) in ((5, 1), (7, 1), (5, 2), (7, 2), (5, 3)):
        K, mx = exact_K(b, s)
        print(f"(b,s)=({b},{s}): exact K = {float(K):10.3f}   max single dual norm = {float(mx):8.3f}   "
              f"recursive bound = {float(K_bound(b, s)):10.3f}   ratio = {float(K_bound(b, s) / K):.2f}")
