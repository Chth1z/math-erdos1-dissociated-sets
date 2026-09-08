# Exact certificate format and reproduction

Run `python experiments/make_certificates.py` to generate all 15 published rows,
including `(17,4)`, and `python tools/verify_certificates.py` to verify the witnesses.
Python 3.12 or later and the repository requirements are sufficient. The independent
verifier does not import any production module from `src/`.

`certificates.json` is a list of `tensor-certificate-v2` rows. The strings `f_bound`
and `N_over_2^n` are the exact six-place upward roundings of the constructed set's
ratios. Fields ending in `_approx` are informational floating-point approximations.
`K` and `Delta_exact` are exact rational strings. The numerator of the exact f ratio
is reconstructed as the product of level maxima; its denominator is
`2 ** f_denominator_exponent`. No floating-point number participates in acceptance
of a mathematical upper bound.

Each row links a gzip JSON witness in `witnesses/`. Its
`tensor-adjugate-seeds-v1` integers are canonical positive lowercase hexadecimal
strings, without leading zeroes or a `0x` prefix. All integers are exactly
reconstructible from stored seeds, with no null fields. This format compresses
the adjugate algebraically: two entries determine
all of its entries, rather than duplicating the complete vector.

For a level with `theta`, `rho`, and `a=2*theta+rho`, the layer matrix satisfies

```
(M psi)_i = a*psi_i - theta*psi_(i+1) - theta*psi_(i+2)   (indices modulo b).
```

The witness supplies exact `psi_0`, `psi_1`, `theta`, `rho`, and `det M`. Starting
from the two seeds, the verifier reconstructs coordinates `b-1,...,1` by
`psi_i = theta*(psi_(i+1)+psi_(i+2))/a`, requiring every division to be exact and
the reconstructed `psi_1` to equal its seed. This recovers the entire vector.

The determinant is computed independently using `T_0=2`, `T_1=-theta`,
`T_j=-theta*T_(j-1)+a*theta*T_(j-2)`, and, for odd b,
`det M=a**b+T_b-theta**b`. This follows by taking the power sums of the two roots
of `theta*x*x+theta*x-a`. The verifier then checks the remaining coordinate-zero
identity. Since the determinant is nonzero, this establishes that the recovered
vector is the actual adjugate column, with its correct scale. It also checks its
gcd, positivity, the interlevel recursion, independently computed K and Delta,
the theorem's Q hypotheses, the exact shape inequality, and exact upward rounding.

Hash encodings are fixed as follows; every SHA-256 is a full 64-character digest:

- `sha256_witness_gzip`: the bytes of the compressed file.
- `sha256_witness_json`: the uncompressed ASCII JSON bytes.
- `sha256_psi`: all recovered psi entries in lowercase hex, each followed by LF.
- `sha256_max_u`: ASCII decimal digits of max u, without a trailing newline.

Hashes check consistency of the published material; mathematical verification
comes from the exact identities and gcd checks, not the hashes alone. Compressed
witnesses can be regenerated with the generator; gzip header differences between
Python versions do not change the mathematical values or decoded-vector digests.

`experiments/exp09_certificate_checks.py` compares the production Cayley-Hamilton
algorithm, the retained Bareiss implementation, direct SymPy Berkowitz adjugates,
and witness reconstruction on 20 deterministic test matrices. It also checks
decimal/logarithm boundary cases, rejects corrupted witnesses, and generates
`first_level_resultants.json`. That file includes the full first-level adjugate
polynomials, determinant polynomials, a1/c1 polynomials, exact nonzero resultants,
and factorizations for b=5,7,11,13. Their polynomial multiplication identities are
verified exactly before the resultants are accepted.

The chosen l controls a target before rounding Q to its progression and searching
for a primitive instance. It does not by itself guarantee a loss below exp(.005).
For every published row the independent verifier checks the stronger rational
statement `1/(1-(d-1)*(Q-2**l)/2**l) < 1.006`, which bounds the actual factor
`(Q/2**l)**(d-1)` and is itself smaller than `exp(.006)`.
