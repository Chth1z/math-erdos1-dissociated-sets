"""Produce the certificate table used in the paper (results/certificates.json, results/summary.md).

For each (b, s): choose l = ceil(log2(200 (d-1) K)), find the smallest good Q = 2^{s+1} t >= 2^l + K,
and record the exact bound  f(l d) = max(u) / 2^{l(d-1)}  together with the level data.
Everything is exact integer arithmetic; the certificate is "content(psi_k) = 1 for all k".
"""
import hashlib
import json
import sys
from fractions import Fraction
from math import log2
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.set_int_max_str_digits(0)

from tensor import K_bound, delta, tensor_certificate  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PARAMS = [(5, 1), (7, 1), (11, 1), (13, 1), (5, 2), (7, 2), (11, 2), (13, 2), (5, 3), (7, 3), (11, 3), (13, 3), (17, 3), (13, 4)]
if "--big" in sys.argv:
    PARAMS.append((17, 4))


def sha(n: int) -> str:
    return hashlib.sha256(str(n).encode()).hexdigest()[:16]


def main():
    rows = []
    for (b, s) in PARAMS:
        d = b ** s
        K = K_bound(b, s)
        Delta = delta(b, s)
        l = int(log2(200 * (d - 1) * float(K))) + 1
        Q0 = 2 ** (s + 1)
        t = int((2 ** l + K) // Q0) + 1
        tries = 0
        while True:
            tries += 1
            cert = tensor_certificate(b, s, Q0 * t)
            if cert.primitive:
                break
            t += 1
        Q = Q0 * t
        umax, umin = cert.u_max(), cert.u_min()
        f_exact = Fraction(umax, 2 ** (l * (d - 1)))
        ratio_Q = Fraction(umax, Q ** (d - 1))
        row = {
            "b": b, "s": s, "d": d, "l": l, "n": l * d, "Q": Q, "tries": tries,
            "K": str(K), "K_float": float(K), "Delta": float(Delta),
            "f_bound": float(f_exact), "maxu_over_Q^(d-1)": float(ratio_Q),
            "N_over_2^n": float(f_exact) / 2,
            "min_u_over_max_u": float(Fraction(umin, umax)),
            "digits_max_u": len(str(umax)), "sha256_max_u": sha(umax),
            "levels": [
                {"k": L.k, "cont": L.cont, "digits_theta": len(str(L.theta)), "digits_rho": len(str(L.rho)),
                 "theta": str(L.theta) if len(str(L.theta)) <= 200 else None,
                 "rho": str(L.rho) if len(str(L.rho)) <= 200 else None,
                 "u": [str(x) for x in L.u] if len(str(max(L.u))) <= 60 else None,
                 "a_digits": len(str(L.a)), "c_digits": len(str(L.c))}
                for L in cert.levels
            ],
        }
        rows.append(row)
        print(f"({b},{s}) d={d} l={l} n={l*d} Q={Q} f<={float(f_exact):.6f} N/2^n<={float(f_exact)/2:.6f} "
              f"Delta={float(Delta):.5f} conts={[L.cont for L in cert.levels]}", flush=True)
    out = ROOT / "results" / "certificates.json"
    out.write_text(json.dumps(rows, indent=1), encoding="utf-8")
    md = ["| b | s | d=b^s | l | n=ld | Q | f(n) bound = max u / 2^{l(d-1)} | N/2^n bound | Delta_{b,s} | K_s |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['b']} | {r['s']} | {r['d']} | {r['l']} | {r['n']} | {r['Q']} | {r['f_bound']:.6f} | "
                  f"{r['N_over_2^n']:.6f} | {r['Delta']:.5f} | {r['K_float']:.1f} |")
    (ROOT / "results" / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("written results/certificates.json and results/summary.md")


if __name__ == "__main__":
    main()
