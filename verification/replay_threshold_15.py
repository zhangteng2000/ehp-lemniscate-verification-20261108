"""Replay all scalar verifications for arbitrary degrees n >= 15."""
from pathlib import Path
import subprocess,sys,json

if __name__=="__main__":
    if not __debug__ or sys.flags.optimize:raise SystemExit("Run with assertions enabled.")
    here=Path(__file__).resolve().parent
    for name in ("effective_line_constants.py","audit_effective_line_256.py","audit_n256_identities.py",
                 "verify_interpolation_degrees.py","independent_interpolation_check.py",
                 "verify_peeled_line.py","small_degree_circle_certificate.py","audit_peeled_line.py"):
        subprocess.run([sys.executable,"-S","-B",str(here/name)],cwd=here,check=True)
    def data(n):return json.loads((here/n).read_text())
    assert [r["n"] for r in data("interpolation_degree_audit.json")]==[r[0] for r in data("independent_interpolation_results.json")]==list(range(28,256))
    assert [r[0] for r in data("peeled_primary.json")["checks"]]==[r["n"] for r in data("small_degree_circle_intervals.json")]==list(range(15,28))
    assert sorted(map(int,data("peeled_independent.json")))==list(range(15,28))
    print("Both finite strips and the unbounded analytic range have passed their constituent verification programs.")

