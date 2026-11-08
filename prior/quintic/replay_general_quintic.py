"""Rebuild and rerun the complete quintic verifier from source in this directory."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import os
import shutil
import subprocess
import sys

here=Path(__file__).resolve().parent
def call(*args):
    subprocess.run([str(s) for s in args],cwd=here,check=True)
def py(*args):call(sys.executable,"-S","-B",*args)
def compile_one(name):
    if os.name=="nt":
        windir=Path(os.environ.get("WINDIR",r"C:\Windows"))
        candidates=[windir/"Microsoft.NET/Framework64/v4.0.30319/csc.exe",shutil.which("csc")]
        compiler=next((str(p) for p in candidates if p and Path(p).is_file()),None)
        if not compiler:raise SystemExit("A C# compiler is required (the 64-bit .NET Framework csc.exe is supported).")
        call(compiler,"/nologo","/checked+","/optimize+","/platform:x64","/out:"+name+".exe",name+".cs")
    else:
        compiler=shutil.which("mcs")
        if not compiler or not shutil.which("mono"):raise SystemExit("Install the Mono C# compiler and runtime (mcs and mono).")
        call(compiler,"-checked+","-optimize+","-platform:x64","-out:"+name+".exe",name+".cs")
def verify_published():
    manifest=here/"published_sha256.json"
    if not manifest.exists():return
    for n,h in json.loads(manifest.read_text()).items():
        actual=sha256((here/n).read_bytes()).hexdigest()
        if h!=actual:raise SystemExit("Published checksum mismatch: "+n)
if __name__=="__main__":
    if not __debug__:raise SystemExit("Run without -O: the verifiers use assertions.")
    parser=argparse.ArgumentParser();parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--second-source",action="store_true",help="also rerun the rectangular traversal with NumPy and unbounded integers")
    args=parser.parse_args();verify_published()
    for name in ("audit_quintic_identities.py","audit_quintic_rational.py","audit_phase_roundoff.py"):py(name)
    compile_one("QuinticPhaseGrid");compile_one("QuinticPhaseAudit")
    for k in "AB":py("quintic_universal_certificate.py",k,"--workers",args.workers,"--replay")
    py("audit_quintic_phase_tree.py","--workers",args.workers)
    py("independent_quintic_shapes.py","--out","quintic_shape_disk_audit.json")
    if args.second_source:
        call(sys.executable,"-B","quintic_shape_integer_cover.py","--unbounded","--log","quintic_shape_arbitrary_audit.json")
        a=json.loads((here/"quintic_shape_arbitrary_audit.json").read_text())
        b=json.loads((here/"quintic_shape_integer_final.json").read_text())
        assert a==b and a["residual_count"]==0
    print("All requested certificates recomputed. Both phase domains and the complete source domain are covered.")

