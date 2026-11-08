"""A second vertex evaluation and an independent recursive proof of coverage."""
from pathlib import Path
from itertools import product
from hashlib import sha256
import argparse
import json
import subprocess
import os

here=Path(__file__).resolve().parent
root=here/"quintic_universal_phase_certificate"
meshes=((8,16,32,256,128),(16,32,64,512,256),(32,64,128,512,256))

def run(kind,workers):
    bound=7 if kind=="A" else 5
    covered=set()
    levels=[]
    def passes(n,t,h):
        # lower area = 333/106 * (5/2 or 3) * n/(t*h).
        return 333*(10 if kind=="A" else 6)*n>106*bound*t*h
    for i,(q,a,b,t,h) in enumerate(meshes):
        source=root/f"{kind}_{i}_counts.txt"
        if not source.exists():break
        points=[]
        for row in source.read_text().splitlines():
            u,v,w,c=map(int,row.split())
            if passes(c,t,h):points.append((u,v,w))
        path=root/f"{kind}_{i}_leaf_input.txt"
        path.write_text("".join(f"{u} {v} {w}\n" for u,v,w in points),encoding="ascii")
        dest=root/f"{kind}_{i}_audit.txt"
        runner=[] if os.name=="nt" else ["mono"]
        subprocess.run([*runner,str(here/"QuinticPhaseAudit.exe"),kind,*map(str,(q,a,b,t,h)),str(path),str(dest),
                        str(here/"quintic_sine_independent.txt"),str(workers)],check=True)
        rows=list(map(lambda z:tuple(map(int,z.split())),dest.read_text().splitlines()))
        assert len(rows)==len(set(points))==len(points)
        assert [r[:3] for r in rows]==sorted(points)
        assert all(passes(v,t,h) and 0<=v<=t*h for *_,v in rows)
        for u,v,w in points:
            assert 0<=u<q and 0<=v<a and 0<=w<b
            covered.add((i,u,v,w))
        levels.append(dict(mesh=[q,a,b,t,h],leaves=len(points),min=min(r[3] for r in rows),
                           sha256=sha256(dest.read_bytes()).hexdigest()))
        print(kind,i,levels[-1],flush=True)
    consumed=0
    def visit(i,u,v,w):
        nonlocal consumed
        if (i,u,v,w) in covered:
            consumed+=1
            return
        assert i+1<len(levels),(i,u,v,w)
        for d,e,f in product((0,1),repeat=3):visit(i+1,2*u+d,2*v+e,2*w+f)
    for u,v,w in product(range(8),range(16),range(32)):visit(0,u,v,w)
    assert consumed==len(covered)
    return levels

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--workers",type=int,default=2);a=p.parse_args()
    out={k:run(k,a.workers) for k in "AB"}
    (root/"independent_phase_audit.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")

