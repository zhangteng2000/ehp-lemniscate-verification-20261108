"""Independent one-coefficient circle comparison for all degrees 9--27.

The constants use a positive binomial sum and a split beta series. The
interval cover is by trisection, fractional powers use integer bisection,
and the curvature bound uses a separate row-tail estimate.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
import argparse,json
from audit_low_degree_transfer import const
from audit_peeled_line import PLO,PHI,root_floor,atan_inv
from audit_capacity_transfer import curvature
from audit_quintic_rational import I

def run(n):
    assert 9<=n<=27
    x,b,Il,Ih,gap,H=const(n)
    L=Q(41,50) if n==9 else Q(7,8) if n<=14 else Q(49,60)
    t=Q(1,320) if n==9 else Q(1,(8 if n<=14 else 128)*n*n)
    loss=None
    if n<=14:
        U=I(2*PLO*(1+(n-1)*Il/4),2*PHI*(1+(n-1)*Ih/4))
        c=(x,b,None,I(Il,Ih),U,None)
        hn=Q(9,20)+Q(4,9)*t if n==9 else None
        goal=Q(4,5) if n==9 else Q(5,8)
        loss=curvature(n,goal,t,c,hn_bound=hn)
    start=2*t/(1+t*t);stack=[(start,Q(1,2)),(Q(1,2),Q(1))];done=[]
    while stack:
        lo,hi=stack.pop();gain=Q(0)
        if hi<=Q(1,2):
            c=root_floor(1-hi*hi,2);K=Q(9,10) if n==9 else H
            gain=2*PLO*(Il-1)*max(0,b*c-K*hi)**2
        if gain+L*root_floor(lo,n)>gap:done.append((lo,hi));continue
        assert hi-lo>Q(1,2**28),(n,str(lo),str(hi))
        a=(2*lo+hi)/3;z=(lo+2*hi)/3
        stack.extend(((lo,a),(a,z),(z,hi)))
    done.sort();assert done[0][0]==start and done[-1][1]==1 and all(a[1]==b[0] for a,b in zip(done,done[1:]))
    print('Independent',n,'one coefficient,',len(done),'closed intervals; curvature',loss,flush=True)
    return {'n':n,'tau':str(t),'rows':1,'R9_upper':('9/10' if n==9 else None),'curvature_upper_display':loss,
            'intervals':[[str(a),str(b)] for a,b in done]}

if __name__=='__main__':
    a,b=atan_inv(5);c,d=atan_inv(239)
    assert PLO<16*a-4*d<16*b-4*c<PHI
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=list(range(9,28)));a=p.parse_args()
    out={str(n):run(n) for n in a.n}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8', newline='\n')
