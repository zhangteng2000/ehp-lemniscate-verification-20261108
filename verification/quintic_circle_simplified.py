"""The quintic circle comparison with three rows and six diameter intervals.

The line estimate is supplied by the separate quintic line certificate.
This file replaces the old 25-row circle comparison; assertions are exact.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as F
import json
from audit_quintic_rational import I,PI,taylor,vin
from coherent_projection_transfer import constants
from diagonal_curvature_bound import certified

def verify():
    x,b,C,moment,U,_=constants(5)
    H=(U/(2*PI)-1)**F(1,2)
    hn=(1/(8*I(2)**F(1,2))).h+F(2,5)*F(1,120)
    curvature=certified(5,F(3),F(1,120),F(2,5),constants(5),hn_bound=hn,components=True)
    assert curvature['total']<F(9,10)<F(17,16)
    normal=[];sums=[]
    for m in range(1,4):
        P=[];value=F(0)
        for k in range(m):
            value+=taylor(-b,k)*taylor(b,m-k);P.append(value)
        assert value+taylor(-b,m)==0
        row=[sum(P[i]*taylor(-b,abs(i-j))/taylor(1+b,abs(i-j))
                 for i in range(m)) for j in range(m)]
        def pair(q):return moment*sum(v*w for v,w in zip(row,q))-P[0]*q[0]
        N=pair(P);assert N.l>0;normal.append(N);band={}
        for J in range(1,m+1):
            for l in range(m-J+1):
                K=(F(1,5)+1/(4*I(2)**F(1,2)) if (J,l)==(1,0) else
                   H*F((J+l)**(J+l),J**J*(l**l if l else 1)))
                band[J+l]=band.get(J+l,F(0))+pair(vin(m,J,l,b)).absupper()*K.h
        sums.append(band)
    ds=[F(480,14401)]+[F(1,k) for k in (24,16,10,6,4,3)]
    truncations=(3,3,3,3,2,1)
    margins=[]
    for lo,hi,M in zip(ds,ds[1:],truncations):
        e=hi/2;head=I(0)
        for m in range(M):
            head+=(1-e*e)**(m+1)*normal[m]-2*sum(v*e**j for j,v in sums[m].items())
        margin=F(17,16)*I(lo/2)**x+2*PI*head-(U-C)
        assert margin.l>F(3,1000)
        margins.append({'left':str(lo),'right':str(hi),'rows':M,'margin_gt':'3/1000'})
    assert (F(17,16)*I(F(1,6))**x-(U-C)).l>F(1,50)
    out={'status':'all exact checks passed','tau':'1/120','R':'3','b0':'2/5',
         'curvature_lt':'9/10','max_rows':3,'diameter_intervals':margins}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Quintic: curvature <9/10; three rows on six closed intervals passed.')
    return out

if __name__=='__main__':assert __debug__;verify()
