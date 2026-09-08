"""Cross-check independent exact algorithms and publish first-level resultants."""
from __future__ import annotations
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import random
import sys
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from certificates import ceil_log2, upper_decimal
from tensor import adjugate_element, conj, level_adjugate

spec = importlib.util.spec_from_file_location("independent_verifier", ROOT / "tools/verify_certificates.py")
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


def arithmetic_checks():
    values = [Fraction(0), Fraction(1, 10**6), Fraction(10**40+1,10**46),
              Fraction(20607218983907535,10**17), Fraction(1,3), Fraction(2,3)]
    for value in values:
        gap = Fraction(upper_decimal(value)) - value
        assert 0 <= gap < Fraction(1, 10**6)
    assert upper_decimal(Fraction(20607218983907535,10**17)) == "0.206073"
    for exponent in range(-100,101):
        power = Fraction(2) ** exponent
        assert ceil_log2(power) == exponent
        assert ceil_log2(power * Fraction(1000001,1000000)) == exponent + 1
        assert ceil_log2(power * Fraction(999999,1000000)) == exponent
    print("PASS exact upper rounding and log2 boundaries (including 201 exact powers of two)", flush=True)
    rng = random.Random(20260908)
    for b in (3,5,7,9,11):
        for _ in range(4):
            theta, rho = rng.randrange(1,50), rng.randrange(1,30)
            psi, determinant = level_adjugate(b, theta, rho)
            phi = [2*theta+rho,-theta,-theta] + [0]*(b-3)
            assert adjugate_element(conj(phi)) == (psi, determinant)
            M = sp.Matrix(b,b,lambda i,j: (2*theta+rho if i == j else 0)
                          - (theta if j == (i+1)%b else 0) - (theta if j == (i+2)%b else 0))
            assert list(M.adjugate(method="berkowitz")[:,0]) == psi
            assert M.det(method="berkowitz") == determinant
            assert verifier.reconstruct_adjugate(b,theta,rho,psi[:2]) == (psi,determinant)
            try:
                verifier.reconstruct_adjugate(b,theta,rho,[psi[0]+1,psi[1]])
            except ValueError:
                pass
            else:
                raise AssertionError("the verifier accepted a corrupted adjugate seed")
    print("PASS 20 cases: Cayley-Hamilton, Bareiss, Berkowitz, and witness reconstruction; corrupt seeds rejected", flush=True)


def resultant_certificates():
    q = sp.Symbol("q")
    result = []
    for b in (5,7,11,13):
        samples = [level_adjugate(b, value, 1)[0] for value in range(1,b+1)]
        psi = [sp.Poly(sp.interpolate([(value,samples[value-1][i]) for value in range(1,b+1)],q),q)
               for i in range(b)]
        # Independently certify the interpolated polynomial adjugate by its full
        # multiplication identity and a scalar determinant recurrence.
        a = sp.Poly(2*q+1,q)
        T0, T1 = sp.Poly(2,q), sp.Poly(-q,q)
        for _ in range(2,b+1):
            T0,T1 = T1,-sp.Poly(q,q)*T1+a*sp.Poly(q,q)*T0
        determinant = a**b + T1 - sp.Poly(q**b,q)
        for i in range(b):
            assert a*psi[i] - sp.Poly(q,q)*(psi[(i+1)%b]+psi[(i+2)%b]) == (determinant if i == 0 else sp.Poly(0,q))
        a1, c1 = psi[0], 2*psi[1]+psi[2]
        resultant = int(sp.resultant(a1,c1))
        assert resultant != 0 and sp.gcd(a1,c1).degree() == 0
        assert (a1 - determinant + sp.Poly(q/b,q)*determinant.diff()).is_zero
        factors = {str(p): int(e) for p,e in sp.factorint(abs(resultant)).items()}
        entry = {"b":b,"parameter":"q=Q/2", "coefficient_order":"highest degree first",
                 "a1_coefficients":[str(x) for x in a1.all_coeffs()],
                 "c1_coefficients":[str(x) for x in c1.all_coeffs()],
                 "psi_coefficients":[[str(x) for x in p.all_coeffs()] for p in psi],
                 "determinant_coefficients":[str(x) for x in determinant.all_coeffs()],
                 "resultant":str(resultant),"absolute_resultant_factorization":factors}
        result.append(entry)
        print(f"PASS b={b}: Res_q(a1,c1)={resultant}; factors={factors}; polynomial adjugate identity exact", flush=True)
    (ROOT / "results/first_level_resultants.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


def manifest_negative_check():
    row = json.loads((ROOT / "results/certificates.json").read_text(encoding="utf-8"))[0]
    verifier.verify_row(row, ROOT / "results")
    # Lower the actual six-place upper bound by exactly one final-place unit.
    units = int(row["f_bound"].replace(".", "")) - 1
    row["f_bound"] = f"{units // 10**6}.{units % 10**6:06d}"
    try:
        verifier.verify_row(row, ROOT / "results")
    except ValueError as error:
        assert "upward rounding" in str(error)
    else:
        raise AssertionError("the verifier accepted a downward-rounded upper bound")
    print("PASS full manifest verification rejects a bound lowered by one final-place unit", flush=True)


if __name__ == "__main__":
    arithmetic_checks()
    resultant_certificates()
    manifest_negative_check()
