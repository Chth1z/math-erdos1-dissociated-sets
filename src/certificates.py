"""Exact reporting and portable witnesses for the published certificates.

The first two entries determine the cyclic adjugate by a backwards order-two
recurrence. Witness integers use lowercase hexadecimal without a 0x prefix.
"""
from __future__ import annotations
from fractions import Fraction
import hashlib


def ceil_fraction(value: Fraction) -> int:
    value = Fraction(value)
    return -(-value.numerator // value.denominator)


def upper_decimal(value: Fraction, digits: int = 6) -> str:
    """Round a nonnegative rational toward +infinity, without floating point."""
    value = Fraction(value)
    if value < 0 or digits < 0:
        raise ValueError("require a nonnegative value and digit count")
    unit = 10 ** digits
    result = ceil_fraction(value * unit)
    return str(result) if digits == 0 else f"{result // unit}.{result % unit:0{digits}d}"


def ceil_log2(value: Fraction) -> int:
    """Least integer e with 2^e >= value, including exact powers of two."""
    value = Fraction(value)
    if value <= 0:
        raise ValueError("log2 requires a positive value")
    p, q = value.numerator, value.denominator
    exponent = p.bit_length() - q.bit_length()
    enough = p <= q << exponent if exponent >= 0 else p << -exponent <= q
    return exponent if enough else exponent + 1


def integer_sha256(value: int) -> str:
    """SHA-256 of the decimal integer in ASCII, without a final newline."""
    return hashlib.sha256(str(value).encode("ascii")).hexdigest()


def vector_sha256(vector) -> str:
    """SHA-256 of lowercase hex entries, each followed by one ASCII LF."""
    digest = hashlib.sha256()
    for value in vector:
        if value < 0:
            raise ValueError("certificate coordinates must be nonnegative")
        digest.update(f"{value:x}\n".encode("ascii"))
    return digest.hexdigest()


def witness_payload(cert) -> dict:
    return {
        "format": "tensor-adjugate-seeds-v1", "b": cert.b, "s": cert.s, "Q": cert.Q,
        "levels": [
            {"k": L.k, "theta_hex": f"{L.theta:x}", "rho_hex": f"{L.rho:x}",
             "determinant_hex": f"{L.N:x}", "psi_seeds_hex": [f"{L.psi[0]:x}", f"{L.psi[1]:x}"],
             "sha256_psi": vector_sha256(L.psi), "content": L.cont}
            for L in cert.levels
        ],
    }
