"""Recheck every parameter leaf, its ancestors, and the integer slice bounds."""
from collections import Counter
from datetime import datetime,timezone
from fractions import Fraction as F
from pathlib import Path
from math import comb
import argparse
import hashlib
import json
import os
import subprocess
import time
import sympy as sp
from chord_refined import geometry,Outside,Split
from chord_geometry import ctx
from independent_controls import dirs8
from slice_data import sqrtlo,sqrtup

HERE=Path(__file__).resolve().parent
END=datetime(2026,11,7,7,45,tzinfo=timezone.utc).timestamp()


def constants():
    low_x=sqrtlo(sqrtlo(sqrtlo(F(2))))
    high_x=sqrtup(sqrtup(sqrtup(F(2))))
    low_y=sqrtlo(sqrtlo(F(1,2))); high_y=sqrtup(sqrtup(F(1,2)))
    u=F(1);l=r=F(0);M=40
    for i in range(M+1):
        l+=F(comb(2*i,i),8**i)/(i+F(1,8))
        r+=u/(2**i*(i+F(1,2)));u*=F(8*i+7,8*(i+1))
    assert low_x*l+low_y*r>F(11069,1000)
    assert high_x*(l+F(1,2**M))+high_y*(r+F(1,2**M))<12
    z=F(1,256)
    for i in range(1,5):
        assert F(1,4)**(i-1)*(i+(4-i)*z)/(1-z)**2<F(26,25)
    assert (7-3*z)/(1-z)**2<F(8,3)**2
    # The bound in z increases on the required interval.
    assert 3-2*7<0
    assert F(3,5)**2>F(1,3)
    assert F(2)-F(9,32)>0 and (F(2)-F(9,32))**2<3
    assert F(1,6)**4+6*F(1,6)**2+8*F(7,4)*F(1,6)<3
    assert 6*F(1,4)**2+8*F(7,4)*F(1,4)**3+3*F(1,4)**4<1
    assert 10**4-F(16,3)*10**3-2*10**2-16*10-F(35,3)>160*F(64,3)
    c,x,y,z,v=sp.symbols("c x y z v")
    G=lambda a:a**4-4*c*a**3/3-2*a**2+4*c*a
    d=sp.expand(G(1)-G(-1));assert d==16*c/3
    assert sp.factor((G(c)-G(-1))/d-(3-c)*(1+c)**3/(16*c))==0
    assert sp.expand(G(z)-G(-1)-(z**2-1)**2+sp.Rational(4,3)*c*(z-2)*(z+1)**2)==0
    cent=sp.Poly(G(z+c/3)-G(-1),z)
    assert cent.nth(3)==0 and cent.nth(2)==-2-2*c**2/3
    assert cent.nth(0)==1+8*c/3+10*c**2/9-c**4/27
    p=sp.symbols("p")
    prod=sp.expand((p-x-y-z)*(p-x+y+z)*(p+x-y+z)*(p+x+y-z))
    assert sp.expand(prod-(p**4-2*(x*x+y*y+z*z)*p*p-8*x*y*z*p
                          +(x*x+y*y+z*z)**2-4*(x*x*y*y+x*x*z*z+y*y*z*z)))==0
    eq=c**4-6*c*c+(16*v-8)*c-3
    q,a=sp.symbols("q a")
    assert sp.simplify(eq.subs({c:-1+a*q,v:a**3})/(-a**3)-(4*q**3-a*q**4+16-16*a*q))==0
    dirs=dirs8();weights=[]
    for (a,b),(c,d) in zip(dirs,dirs[1:]+[(-dirs[0][0],-dirs[0][1])]):
        assert a*a+b*b==1 and a*d-b*c>0
        weights.append((a*d-b*c)/(1+a*c+b*d))
    return [(w+weights[i-1])/2 for i,w in enumerate(weights)]


class Observer:
    def __init__(self,n,depth,workers,mode,table_path):
        self.n=n;self.uid=0
        binary="CheckSlices" if mode=="independent" else "SliceCounter"
        executable=HERE/(binary+".exe")
        command=[str(executable)] if os.name=="nt" else ["mono",str(executable)]
        self.p=subprocess.Popen(command+[str(table_path),str(depth),str(workers)],
                                stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        assert self.p.stdout.readline().strip()=="READY"
    def send(self,s):
        self.p.stdin.write(s+"\n");self.p.stdin.flush()
        ans=self.p.stdout.readline()
        if not ans:raise RuntimeError("count backend stopped")
        return ans.strip()
    def frame(self,par,data):
        full,target,gens,rad,cm,cr,_,_,_=data
        self.uid+=1
        controls=[self.uid,par,int(full),len(gens)]
        for vec in gens:controls.extend(vec)
        controls+=rad+cm+cr
        a,_=self.send(" ".join(map(str,controls))).split()
        return self.uid,[int(t) for t in a.split(",")]
    def drop(self,i):
        if i>=0:assert self.send("D "+str(i))=="OK"
    def close(self):
        self.p.stdin.close();assert self.p.wait()==0


def expected_tiles():
    y=list(map(F,[0]))+[F(1,4096),F(1,1024),F(1,512),F(1,128),F(1,32),F(1,8),F(9,32)]
    out={}
    for j in range(7):
        assert y[j]<y[j+1]
        for s in range(4):
            for i in range(8):
                out[f"c{j}_h{s}_s{i}"]=(s,(F(2*i+1,16),(y[j]+y[j+1])/2,F(1,16),(y[j+1]-y[j])/2))
    return out


def make_trie(records):
    root={}
    for row in records:
        node=root
        for digit in row[0]:
            assert digit in '0123' and "leaf" not in node
            node=node.setdefault(digit,{})
        assert not node and row[1] in ('P','O')
        node["leaf"]=row
    return root


def audit(certificate,mode,workers,out,group=None,live_deadline=False):
    if mode=="independent":ctx.prec=150
    else:ctx.prec=104
    weights=constants();record=json.loads(certificate.read_text(encoding="utf-8"))
    times=time.monotonic();counts=Counter();margins=[];switches=0;diff=0
    exp=expected_tiles()
    assert set(record["tiles"])==set(exp)
    for n,depth in [(1024,6),(2048,7)]:
        nm=("independent_"+str(n)+".bin" if mode=="independent" else "chord_"+str(n)+"_8.bin")
        obs=Observer(n,depth,workers,mode,HERE/nm)
        for name,(sh,box) in sorted(exp.items()):
            doc=record["tiles"][name]
            assert sh==doc["sheet"] and tuple(map(F,doc["box"]))==box
            assert doc["directions"]==8 and doc["scheme"]=="chord16"
            assert (doc["source_nx"],doc["depth"]) in [(1024,6),(2048,7)]
            if n!=doc["source_nx"] or group is not None and not name.startswith(group):continue
            if live_deadline and time.time()>END:raise TimeoutError("research window")
            tree=make_trie(doc["records"])
            def visit(t,B,parent):
                nonlocal switches,diff
                x,y,a,b=B
                row=t.get("leaf")
                if row is None:
                    assert set(t) in ({"0","1"},{"2","3"})
                ident=-1
                try:
                    res=geometry(sh,x,y,a,b)
                    ident,sums=obs.frame(parent,res)
                    if row is not None:
                        assert row[1]=="P"
                        full,target=res[:2];budget=26 if full else 30
                        upper=F(8,n)*sum(w*(v+budget) for w,v in zip(weights,sums))
                        gap=target-upper
                        if gap<=0:raise ArithmeticError((name,row[0],"nonpositive independent margin",str(gap)))
                        if mode=="primary":assert int(full)==row[2] and sums==row[3] and gap==F(row[4])
                        elif int(full)!=row[2]:switches+=1
                        elif sums!=row[3]:diff+=1
                        counts["positive"]+=1;counts["full" if full else "clipped"]+=1
                        margins.append(gap)
                except Outside:
                    if row is not None:
                        assert row[1]=="O";counts["outside"]+=1
                except (ValueError,ZeroDivisionError,Split):
                    assert row is None
                if row is None:
                    counts["internal"]+=1
                    if "0" in t:
                        visit(t["0"],(x-a/2,y,a/2,b),ident)
                        visit(t["1"],(x+a/2,y,a/2,b),ident)
                    else:
                        visit(t["2"],(x,y-b/2,a,b/2),ident)
                        visit(t["3"],(x,y+b/2,a,b/2),ident)
                obs.drop(ident)
            visit(tree,box,-1);counts["tiles"]+=1
        obs.close()
    report={"status":"PASS","mode":mode,"precision_bits":ctx.prec,"nodes":counts,
            "min_table_coordinate_margin":str(min(margins)),"different_direction_sums":diff,
            "different_target_modes":switches,"elapsed_seconds":time.monotonic()-times,
            "certificate_sha256":hashlib.sha256(certificate.read_bytes()).hexdigest()}
    if group is None:assert counts["tiles"]==224
    if out:out.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return report


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--mode",choices=["primary","independent"],default="independent")
    p.add_argument("--certificate",type=Path,default=HERE/"quartic_atlas_certificate.json")
    p.add_argument("--out",type=Path);p.add_argument("--workers",type=int,default=4)
    p.add_argument("--group");p.add_argument("--live-deadline",action="store_true")
    a=p.parse_args();print(json.dumps(audit(a.certificate,a.mode,a.workers,a.out,a.group,a.live_deadline),indent=2))

