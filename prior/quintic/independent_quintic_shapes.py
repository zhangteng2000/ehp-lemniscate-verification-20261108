"""Independent all-source proof with Gaussian disks and unbounded Python integers."""
from math import isqrt
from collections import Counter
from pathlib import Path
import argparse
import json
S=1<<38


def ceiling(a,d):return -((-a)//d)
def ceilroot(a):
    t=isqrt(a);return t+(t*t!=a)

class Disk:
    __slots__=("a","b","r","upper","lower")
    def __init__(self,a,b=0,r=0):
        self.a,self.b,self.r=a,b,r
        t=isqrt(a*a+b*b)
        self.upper=t+(t*t!=a*a+b*b)+r
        self.lower=max(0,t-r)
    def __add__(self,z):return Disk(self.a+z.a,self.b+z.b,self.r+z.r)
    def __sub__(self,z):return Disk(self.a-z.a,self.b-z.b,self.r+z.r)
    def __mul__(self,z):
        a=self.a*z.a-self.b*z.b;b=self.a*z.b+self.b*z.a
        # |center_A| r_B + |center_B| r_A + r_A r_B, then floor
        # each center component and enclose its two rounding errors.
        e=(self.upper-self.r)*z.r+z.upper*self.r
        return Disk(a//S,b//S,ceiling(e,S)+int(a%S!=0)+int(b%S!=0))
    def rat(self,p,q):
        a=self.a*p;b=self.b*p
        return Disk(a//q,b//q,ceiling(abs(p)*self.r,q)+int(a%q!=0)+int(b%q!=0))
    def conj(self):return Disk(self.a,-self.b,self.r)

one=Disk(S);zero=Disk(0)


def visit(box):
    lr,hr,li,hi,ls,hs,lt,ht=box
    if lr>hs or lr+ls>S:return 1
    def minabs(lo,high):return max(0,lo,-high)
    for a,b,c,d in ((lr,hr,li,hi),(lr-S,hr-S,li,hi),
                     (ls,hs,lt,ht),(ls-S,hs-S,lt,ht),(lr-hs,hr-ls,li-ht,hi-lt)):
        if minabs(a,b)**2+minabs(c,d)**2>S**2:return 1
    def disk(a,b,c,d):
        x=(a+b)//2;y=(c+d)//2
        return Disk(x,y,ceilroot(max(x-a,b-x)**2+max(y-c,d-y)**2))
    u=disk(lr,hr,li,hi);v=disk(ls,hs,lt,ht)
    mu=(one+u+v).rat(1,4)
    pts=[z-mu for z in (zero,one,u,v)]
    squares=[z*z for z in pts]
    f=sum(squares,zero).rat(-1,4)
    g=sum((z*z2 for z,z2 in zip(pts,squares)),zero).rat(-1,12)
    uv=u*v
    # Evaluate the primitive at all four marked critical points directly,
    # using 12F(u)=u^3[-3u^2+5(1+v)u-10v], and its swapped form.
    def at(a,b):
        a2=a*a
        return (a2*a*(a2.rat(-3,1)+(a+uv).rat(5,1)-b.rat(10,1))).rat(1,12)
    vals=[zero,(one.rat(-3,1)+(u+v).rat(5,1)-uv.rat(10,1)).rat(1,12),at(u,v),at(v,u)]
    chords=[a-b for k,a in enumerate(vals) for b in vals[:k]]
    D=min(ceiling(5*S,6),max(q.upper for q in chords))
    if 4*(16*7)**10*f.lower**5 > (17*4)**10*D**2*S**3:return 2
    if 8*(16*5)**15*g.lower**5 > (17*2)**15*D**3*S**2:return 3
    q=max(chords,key=lambda z:z.lower)
    for z in chords:
        det=q.conj()*z
        if abs(det.b)>det.r:return 1
    return 0


def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path);p.add_argument("--max",type=int,default=10000000)
    args=p.parse_args();root=(0,S,0,S,0,S,-S,S);stack=[(root,0)]
    counts=Counter();deep=0
    while stack and sum(counts.values())<args.max:
        box,i=stack.pop();k=visit(box);counts[k]+=1;deep=max(deep,i)
        if k==0:
            # Longest physical coordinate first; reverse tie order relative
            # to the 16-child breadth-first fixed-point implementation.
            _,a=max((box[j+1]-box[j],j) for j in range(0,8,2))
            mid=(box[a]+box[a+1])//2;assert box[a]<mid<box[a+1]
            x=list(box);y=list(box);x[a+1]=mid;y[a]=mid
            stack.extend([(tuple(x),i+1),(tuple(y),i+1)])
    out={"unit":S,"counts":dict(counts),"max_depth":deep,"unprocessed":len(stack)}
    print(out,flush=True)
    if args.out:args.out.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    if stack:raise SystemExit("Incomplete")
if __name__=="__main__":main()
