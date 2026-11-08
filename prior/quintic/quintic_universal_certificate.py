"""Check every parameter box of the universal cubic inequality and its subdivision tree."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import subprocess
import os

HERE=Path(__file__).resolve().parent
ROOT=HERE/"quintic_universal_phase_certificate"
STAGES=[(8,16,32,256,128),(16,32,64,512,256),(32,64,128,512,256),
        (64,128,256,1024,512),(128,256,512,1024,512),(256,512,1024,2048,512)]


def digest(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(1<<18),b""):h.update(block)
    return h.hexdigest()


def levels(kind,limit,workers,replay):
    ROOT.mkdir(exist_ok=True)
    old=None
    stats=[]
    for idx,(q,a,b,K,M) in enumerate(STAGES[:limit]):
        if old is None:
            expected=list(itertools.product(range(q),range(a),range(b))); input_name="all"
        else:
            if not old:break
            expected=sorted(tuple(2*v+d for v,d in zip(p,mask)) for p in old
                            for mask in itertools.product(range(2),repeat=3))
            source=ROOT/f"{kind}_{idx}_input.txt"
            input_name=str(source)
            source.write_text("".join(f"{i} {j} {k}\n" for i,j,k in expected),encoding="ascii")
        target=ROOT/f"{kind}_{idx}_counts.txt"
        if replay or not target.is_file():
            runner=[] if os.name=="nt" else ["mono"]
            subprocess.run([*runner,str(HERE/"QuinticPhaseGrid.exe"),kind,str(q),str(a),str(b),str(K),str(M),
                            input_name,str(target),str(HERE/"quintic_sine_dyadic.txt"),str(workers)],check=True)
        records=[tuple(map(int,s.split())) for s in target.read_text(encoding="ascii").splitlines()]
        assert len(records)==len(expected)
        old=[];mins=None
        for k,(i,j,l,v) in zip(expected,records,strict=True):
            assert (i,j,l)==k and 0<=v<=K*M
            ok=(3330*v>742*K*M if kind=="A" else 1998*v>530*K*M)
            if not ok:old.append(k)
            else:mins=v if mins is None else min(mins,v)
        row={"stage":idx,"mesh":[q,a,b,K,M],"input_count":len(expected),"unresolved":len(old),
             "least_passing_cell_count":mins,"counts_sha256":digest(target)}
        print(kind,row,flush=True);stats.append(row)
    result={"kind":kind,"complete":old==[],"stages":stats,
            "source_sha256":digest(HERE/"QuinticPhaseGrid.cs"),"trig_sha256":digest(HERE/"quintic_sine_dyadic.txt")}
    (ROOT/f"{kind}_summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("kind",choices=["A","B"])
    p.add_argument("--stages",type=int,default=4);p.add_argument("--workers",type=int,default=2)
    p.add_argument("--replay",action="store_true")
    a=p.parse_args();levels(a.kind,a.stages,a.workers,a.replay)
