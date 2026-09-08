"""Exact checks supporting universal first-level coprimality.

Checks every odd b from 3 through 101, including 3|b controls. This finite
experiment is a regression check of the general proof, not a substitute for it.
No production source module or existing finite-set certificate is modified.

Run: python experiments/exp10_coprimality.py [--max-b 101]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def fibonacci_polynomials(z, limit):
    F = [sp.Poly(0, z), sp.Poly(1, z)]
    for _ in range(2, limit + 1):
        F.append(F[-1] + sp.Poly(z, z) * F[-2])
    return F


def adjugate_column(z, b, F):
    """Compute the complete normalized adjugate by matrix Cayley-Hamilton.

    H v has coordinates v_(i+1)+v_(i+2). For odd b its characteristic
    polynomial is D_b(z)=z^b-L_b(z)-1, L_b=F_(b+1)+z*F_(b-1).
    """
    D = sp.Poly(z ** b - 1, z) - F[b + 1] - sp.Poly(z, z) * F[b - 1]
    coefficients = [int(c) for c in D.all_coeffs()]
    vector = [1] + [0] * (b - 1)
    columns = [[] for _ in range(b)]
    for k in range(b):
        for i, value in enumerate(vector):
            columns[i].append(value)
        vector = [vector[(i+1) % b] + vector[(i+2) % b] for i in range(b)]
        vector[0] += coefficients[k + 1]
    psi = [sp.Poly.from_list(c, z) for c in columns]
    for i in range(b):
        residual = sp.Poly(z, z) * psi[i] - psi[(i+1) % b] - psi[(i+2) % b]
        assert residual == (D if i == 0 else sp.Poly(0, z))
    return psi, D


def monic_coefficients(poly):
    return [str(c) for c in poly.monic().all_coeffs()]


def check_initial_value_formulas(z, t, F, limit):
    """Certify the root parametrization without division by a symbolic root gap."""
    lam = sp.Poly((t*t - 1) / 4, t)
    r, s = sp.Poly((t-1)/2,t), sp.Poly((-t-1)/2,t)
    Y = [None, sp.Poly(1,t), sp.Poly(-2,t)]
    for i in range(3, limit + 1):
        Y.append(lam * Y[i-2] - Y[i-1])
    for i in range(1, limit + 1):
        # t=r-s; multiplication by t also covers the limit t=0 polynomially.
        numerator = (r - 1) * r ** (i-1) - (s - 1) * s ** (i-1)
        assert (sp.Poly(t,t) * Y[i] - numerator).is_zero
        symbolic_y = (-1) ** (i-1) * (F[i] + F[i-1])
        assert (symbolic_y.as_expr().subs(z,sp.Rational(-1,4))
                == (3*i - 2) * sp.Rational(-1,2) ** (i-1))
        assert Y[i].eval(0) == (3*i - 2) * sp.Rational(-1,2) ** (i-1)
        # Check the Fibonacci expression for the same sequence at arbitrary lambda.
        converted = sp.Poly(symbolic_y.as_expr().subs(z,lam.as_expr()),t)
        assert (converted - Y[i]).is_zero
    return Y, r, s


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-b", type=int, default=101)
    args = parser.parse_args()
    if args.max_b < 3:
        parser.error("--max-b must be at least 3")
    z, q, t = sp.symbols("z q t")
    start = time.monotonic()
    F = fibonacci_polynomials(z, args.max_b + 2)
    Y, r, s = check_initial_value_formulas(z,t,F,args.max_b + 1)
    rows, messages = [], []
    for b in range(3, args.max_b + 1, 2):
        psi, D = adjugate_column(z,b,F)
        A = sp.Poly(z ** (b-1),z) - F[b-1]
        C = 1 + 2 * F[b] + sp.Poly(z,z) * F[b-1]
        assert psi[0] == A
        assert psi[1] == F[b]
        assert psi[2] == 1 + sp.Poly(z,z) * F[b-1]
        assert C == 2 * psi[1] + psi[2]
        assert (D.diff() - b * A).is_zero
        assert (F[b]**2 - F[b-1]*F[b+1] - sp.Poly(z**(b-1),z)).is_zero
        gcd_z = sp.gcd(A,C).monic()
        expected_z = sp.Poly(z+1,z) if b % 3 == 0 else sp.Poly(1,z)
        assert (gcd_z - expected_z).is_zero
        # Homogenize, using z=2+1/q; this avoids rational-expression simplification.
        Aq = sp.Poly(0,q)
        Cq = sp.Poly(0,q)
        for (degree,), coefficient in A.terms():
            Aq += coefficient * sp.Poly(q ** (b-1-degree) * (2*q+1) ** degree,q)
        for (degree,), coefficient in C.terms():
            Cq += coefficient * sp.Poly(q ** (b-1-degree) * (2*q+1) ** degree,q)
        gcd_q = sp.gcd(Aq,Cq).monic()
        expected_q = sp.Poly(q+sp.Rational(1,3),q) if b % 3 == 0 else sp.Poly(1,q)
        assert (gcd_q - expected_q).is_zero
        assert Aq.eval(0) == 1
        assert A.eval(2) == (2**b + 1) // 3
        # Exact identities giving equivalence of the two periodic conditions.
        assert (Y[b+1] - s*Y[b] - (r-1)*r**(b-1)).is_zero
        assert (Y[b+1] - r*Y[b] - (s-1)*s**(b-1)).is_zero
        assert Y[b].eval(0) != 0
        row = {"b":b,"gcd_z_monic":str(gcd_z.as_expr()),
               "gcd_q_monic":str(gcd_q.as_expr()),"coprime":b % 3 != 0,
               "degree_A_z":A.degree(),"degree_C_z":C.degree(),
               "A_z_coefficients":[str(c) for c in A.all_coeffs()],
               "C_z_coefficients":[str(c) for c in C.all_coeffs()],
               "A_q_coefficients":[str(c) for c in Aq.all_coeffs()],
               "C_q_coefficients":[str(c) for c in Cq.all_coeffs()],
               "checks":{"complete_adjugate_identity":True,"three_entry_formula":True,
                         "derivative_formula":True,"cassini":True,"root_periodicity_equivalence":True,
                         "repeated_root_excluded":True}}
        rows.append(row)
        message = f"PASS b={b}: gcd_z={gcd_z.as_expr()}, gcd_q={gcd_q.as_expr()}; adjugate, Cassini, periodicity, repeated root"
        messages.append(message)
        print(message,flush=True)
    result = {"python":platform.python_version(),"sympy":sp.__version__,
              "coefficient_order":"highest degree first","normalization":"z=(2q+1)/q; A_q=q^(b-1)*A_z(z), likewise C",
              "range":{"min_b":3,"max_b":args.max_b,"odd_cases":len(rows)},
              "general_formula":{"F":"F_0=0, F_1=1, F_(n+1)=F_n+z*F_(n-1)",
                                 "A_z":"z^(b-1)-F_(b-1)","C_z":"1+2*F_b+z*F_(b-1)",
                                 "y_i":"(-1)^(i-1)*(F_i+F_(i-1))",
                                 "repeated_root":"lambda=-1/4: y_i=(3i-2)*(-1/2)^(i-1)"},
              "initial_value_formula_checks":args.max_b+1,
              "elapsed_seconds":round(time.monotonic()-start,3),"rows":rows}
    destination = ROOT / "results/coprimality_checks.json"
    destination.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    final = (f"ALL PASS: {len(rows)} odd b; gcd=1 exactly when 3 does not divide b; "
             f"Binet/repeated-root checks for i=1..{args.max_b+1}; {result['elapsed_seconds']:.3f}s. "
             "Finite checks support the separate general proof.")
    messages.append(final)
    (ROOT / "results/coprimality_checks.log").write_text("\n".join(messages)+"\n",encoding="utf-8")
    print(final,flush=True)


if __name__ == "__main__":
    main()
