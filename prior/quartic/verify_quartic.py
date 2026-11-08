"""Build both source tables and both integer counters, then verify the full cover."""
from argparse import ArgumentParser
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
EXPECTED = "964241c49502d2d7c26dad063cca2edf679fb1d0ac64a61379fe0fc70b1059d9"
TABLE_HASH = {
    1024: "fbfe0c110f774d222eff81d670835892bc9c45d368efe93dea4d5be59be832c4",
    2048: "2073a8d86fa9099db87cc5a8ae279eb3e210f2ad17631d6c2393f3c6b5434fe8",
}
MARGIN = {
    "primary": Fraction(54996863, 7902739824640),
    "independent": Fraction(899883551, 7902739824640),
}
def run(args):
    subprocess.run([str(arg) for arg in args], cwd=HERE, check=True)

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def compile_counters():
    if os.name == "nt":
        default = Path(os.environ["SystemRoot"]) / "Microsoft.NET/Framework64/v4.0.30319/csc.exe"
        compiler = os.environ.get("CSC") or (str(default) if default.exists() else shutil.which("csc"))
    else:
        compiler = os.environ.get("CSC") or shutil.which("mcs")
        if not shutil.which("mono"):
            raise RuntimeError("Install the 64-bit Mono runtime and C# compiler.")
    if not compiler:
        raise RuntimeError("No C# compiler found. Set CSC to its absolute path.")
    for name in ("SliceCounter", "CheckSlices"):
        if os.name == "nt":
            run([compiler, "/nologo", "/target:exe", "/optimize+", "/checked+", "/platform:x64",
                 "/out:" + name + ".exe", name + ".cs"])
        else:
            run([compiler, "-optimize+", "-checked+", "-platform:x64",
                 "-out:" + name + ".exe", name + ".cs"])

def verify(workers):
    if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
        raise RuntimeError("Use an unoptimized Python interpreter. Assertions are proof checks.")
    if sys.maxsize <= 2**32 or workers < 1:
        raise ValueError("A 64-bit interpreter and a positive worker count are required.")
    cert = HERE / "quartic_atlas_certificate.json"
    if digest(cert) != EXPECTED:
        raise ValueError("Certificate hash mismatch.")
    compile_counters()
    for n, expected in TABLE_HASH.items():
        first = HERE / f"chord_{n}_8.bin"
        second = HERE / f"independent_{n}.bin"
        run([sys.executable, "slice_data.py", first, "--nx", n, "--directions", 8, "--family", "chord"])
        run([sys.executable, "independent_controls.py", second, n])
        if digest(first) != expected or digest(second) != expected:
            raise ValueError(f"Table mismatch for {n} midpoints.")
    for mode in MARGIN:
        report = HERE / f"replayed_{mode}.json"
        run([sys.executable, "audit_atlas.py", "--mode", mode, "--workers", workers, "--out", report])
        data = json.loads(report.read_text(encoding="utf-8"))
        counts = data["nodes"]
        if (data["status"] != "PASS" or data["certificate_sha256"] != EXPECTED or
            counts["tiles"] != 224 or counts["outside"] != 646 or counts["positive"] != 20809 or
            counts["clipped"] != 4211 or counts["full"] != 16598 or
            Fraction(data["min_table_coordinate_margin"]) < MARGIN[mode]):
            raise ArithmeticError(f"Exhaustive verification failed in {mode} mode.")
    print("PASS: complete quartic cover; both integer kernels, both source tables and both margins.")
if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=3)
    verify(parser.parse_args().workers)

