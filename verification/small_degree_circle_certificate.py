"""Dyadic-interval proof of the scalar circle comparison for n=15,...,27."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as Q
from audit_quintic_rational import I,PI,J
import json
def check(n):
    x=Q(1,n);b=(1-x)/2;L=Q(49,60)
    I0=I(2)**(2-x)*J(2-x)/PI
    Hsq=(n-1)*I0/4;H=(Hsq**Q(1,2)).h
    U=2*PI*(1+Hsq)
    C=I(2)**(1+x)*J(x)
    gap=(U-C).h
    v=Q(1);ten=Q(0)
    for k in range(1,11):
        v=v*(b-k+1)/k;ten+=v*v
    low=Q(1,128*n*n);todo=[(low,Q(1,2)),(Q(1,2),Q(1))];leaves=[]
    while todo:
        lo,hi=todo.pop();extra=Q(0)
        if hi<=Q(1,2):
            c=(I(1-hi*hi)**Q(1,2)).l
            extra=2*PI.l*ten*max(Q(0),b*c-H*hi)**2
        need=(gap-extra)/L
        if need<=0 or need**n<lo:
            leaves.append([str(lo),str(hi)])
        else:
            assert hi-lo>Q(1,2**30),("uncovered",n,float(lo),float(hi))
            mid=(lo+hi)/2;todo.extend([(mid,hi),(lo,mid)])
    leaves.sort(key=lambda u:Q(u[0]))
    assert Q(leaves[0][0])==low and Q(leaves[-1][1])==1
    assert all(a[1]==b[0] for a,b in zip(leaves,leaves[1:]))
    return {"n":n,"intervals":leaves,"g_upper_display":float(gap)}
if __name__=="__main__":
    assert __debug__
    out=[check(n) for n in range(15,28)]
    print("Closed every real chord diameter for every degree 15 through 27",[(r["n"],len(r["intervals"])) for r in out])
    Path(__file__).with_name("small_degree_circle_intervals.json").write_text(json.dumps(out,indent=1)+"\n")

