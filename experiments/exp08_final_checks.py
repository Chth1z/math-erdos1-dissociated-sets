"""Experiment 8: end-to-end brute-force confirmation for s=2, true dissociation threshold
for small cases, and certificates for larger (b,s)."""
import sys
import json
import importlib.util
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from dissociated import is_Q_dissociated, is_dissociated_big as is_dissociated, lift, f_ratio  # noqa: E402
from tensor import K_bound, tensor_certificate  # noqa: E402

sys.set_int_max_str_digits(0)


def brute_force_s2():
    print("== brute force: (b,s)=(5,2), l=1, |A|=25 ==")
    b, s = 5, 2
    K = K_bound(b, s)
    Q = 2 ** s
    found = 0
    while found < 3:
        Q += 2 ** (s + 1)
        if Q < 2 + K:
            continue
        cert = tensor_certificate(b, s, Q)
        if cert.primitive:
            u = cert.full_u()
            A = lift(u, 1)
            print(f"Q={Q}: primitive, |A|={len(A)}, dissociated={is_dissociated(A)}, f(A)={f_ratio(A):.4f}, "
                  f"min u={min(u)} max u={max(u)}")
            found += 1


def true_threshold():
    print("\n== true dissociation threshold vs guaranteed Q-K (b=5, s=1) ==")
    b = 5
    K = K_bound(b, 1)
    for Q in (16, 24, 40):
        cert = tensor_certificate(b, 1, Q)
        if not cert.primitive:
            continue
        u = cert.full_u()
        # Monotonicity in the relation radius makes these two checks sufficient;
        # do not repeat the whole enumeration for every smaller radius.
        assert is_Q_dissociated(u, Q)
        assert not is_Q_dissociated(u, Q + 1)
        print(f"Q={Q}: u={u}, guaranteed (Q-K)={float(Q - K):.2f}, true threshold Q'={Q}")


def big_certificates():
    print("\n== independent exact witness verification for larger (b,s) ==")
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location("independent_verifier", root / "tools/verify_certificates.py")
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    rows = json.loads((root / "results/certificates.json").read_text(encoding="utf-8"))
    for (b, s) in ((17, 3), (13, 4), (17, 4)):
        row = next(r for r in rows if (r["b"],r["s"]) == (b,s))
        verifier.verify_row(row, root / "results")
        print(f"(b,s)=({b},{s}) d={row['d']} l={row['l']} n={row['n']}: Q={row['Q']} "
              f"(recorded tries={row['tries']}) f(n) <= {row['f_bound']} "
              f"N/2^n <= {row['N_over_2^n']} conts={[L['cont'] for L in row['levels']]}")


if __name__ == "__main__":
    brute_force_s2()
    true_threshold()
    big_certificates()
