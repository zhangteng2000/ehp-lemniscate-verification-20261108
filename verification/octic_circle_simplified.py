"""The degree-eight circle comparison with only two exterior rows."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
import json
from coherent_projection_transfer import constants,rows
from diagonal_curvature_bound import certified
from audit_quintic_rational import I,PI
from audit_octic_majorant_exact import verify as majorant

def verify():
    majorant()
    n=8;L=Q(39,50);K=Q(53,100);t=Q(1,300);M=2
    cc=constants(n);x,b,C,I0,U,J=cc;gap=(U-C).h
    h=K/2+Q(7,16)*t
    curve=certified(n,3,t,Q(2,5),cc,hn_bound=h)
    assert curve<Q(7,10)<L
    data=rows(n,M,cc);start=2*t/(1+t*t)
    stack=[(start,Q(3,4)),(Q(3,4),Q(1))];done=[]
    while stack:
        lo,hi=stack.pop();gain=Q(0)
        if hi<=Q(3,4):
            c=(I(1-hi*hi)**Q(1,2)).l;cb=(I(1-lo*lo)**Q(1,2)).h
            for m,(N,band) in enumerate(data,1):
                error=sum(v*hi**j*cb**(m-j) for j,v in band.items())
                value=max(Q(0),c**m*N.l-error)**2/N.h
                if m==1:value=max(value,(I0.l-1)*max(Q(0),b*c-K*hi)**2)
                gain+=value
        # Pure Jensen is enough: no angular-response calculation is used here.
        target=(gap-2*PI.l*gain)/L
        if target<=0 or target**n<lo:done.append((lo,hi));continue
        assert hi-lo>Q(1,2**24),('diameter not closed',str(lo),str(hi))
        mid=(lo+hi)/2;stack.extend(((lo,mid),(mid,hi)))
    done.sort()
    assert done[0][0]==start and done[-1][1]==1 and all(a[1]==b[0] for a,b in zip(done,done[1:]))
    out={'status':'all exact checks passed','L':str(L),'K':str(K),'tau':'1/300',
         'R':'3','b0':'2/5','curvature_lt':'7/10','rows':M,'angular_lower':'Jensen',
         'intervals':[[str(a),str(b)] for a,b in done]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Degree eight:',M,'rows,',len(done),'closed intervals; curvature',float(curve))
    return out

if __name__=='__main__':assert __debug__;verify()
