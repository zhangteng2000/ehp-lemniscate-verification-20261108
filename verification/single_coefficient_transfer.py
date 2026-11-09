"""One inverse coefficient suffices for the circle comparison in degrees 9--27.

All fractional powers and constants are enclosed by outward rational
intervals. The sum over all source-variation coefficients is evaluated
analytically as the single moment I0-1.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
import json,argparse
from audit_quintic_rational import I,PI
from coherent_projection_transfer import constants
from diagonal_curvature_bound import certified

def run(n):
    assert 9<=n<=27
    cc=constants(n);x,b,C,I0,U,_=cc
    L=Q(41,50) if n==9 else Q(7,8) if n<=14 else Q(49,60)
    t=Q(1,320) if n==9 else Q(1,(8 if n<=14 else 128)*n*n)
    H=((U/(2*PI)-1)**Q(1,2)).h
    if n<=14:
        hn=Q(9,20)+Q(4,9)*t if n==9 else None
        loss=certified(n,3,t,Q(2,5),cc,hn_bound=hn)
        assert loss<(Q(4,5) if n==9 else Q(5,8))<L
    else:loss=None  # The analytic curvature256 lemma applies here.
    start=2*t/(1+t*t);stack=[(start,Q(1,2)),(Q(1,2),Q(1))];done=[];least=None
    while stack:
        lo,hi=stack.pop();gain=Q(0)
        if hi<=Q(1,2):
            c=(I(1-hi*hi)**Q(1,2)).l
            K=Q(9,10) if n==9 else H
            gain=2*PI.l*(I0.l-1)*max(0,b*c-K*hi)**2
        target=((U-C).h-gain)/L
        if target<=0 or target**n<lo:
            margin=gain+L*(I(lo)**x).l-(U-C).h
            assert margin>0
            least=margin if least is None else min(least,margin)
            done.append((lo,hi));continue
        assert hi-lo>Q(1,2**28),(n,str(lo),str(hi))
        mid=(lo+hi)/2;stack.extend(((lo,mid),(mid,hi)))
    done.sort();assert done[0][0]==start and done[-1][1]==1 and all(a[1]==b[0] for a,b in zip(done,done[1:]))
    out={'n':n,'tau':str(t),'line_constant':str(L),'rows':1,'R9_upper':('9/10' if n==9 else None),
         'curvature_upper_display':None if loss is None else float(loss),
         'minimum_margin_lower':str(Q((least*(1<<80)).__floor__(),1<<80)),
         'intervals':[[str(a),str(b)] for a,b in done]}
    print(n,'single row,',len(done),'closed intervals; least margin',float(least),flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=list(range(9,28)));a=p.parse_args()
    out={str(n):run(n) for n in a.n}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8', newline='\n')
