"""Exact angular strips where the first derivative of a nonic real slice has a nonreal pair."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial
import json
from pathlib import Path
from verify_conditional_line import up,down
from verify_interpolation_degrees import LOWPI
from verify_peeled_line import variance
G=Path(__file__).resolve().parent
N=9;M=8;U=Q(12,25);L=Q(41,50)
# Upper bounds on the coefficients of z^2,...,z^7 of the centered original F.
CAP={2:Q(103,1000),3:Q(14,1000),4:Q(84,1000),5:Q(311,1000),6:Q(521,1000),7:Q(785,1000)}
def sinlower(a,b):
    if a.__ceil__()<=b:return Q(0)
    d=min(a-a.__floor__(),b.__ceil__()-b)
    t=LOWPI*d
    return max(Q(0),down(sum((-1)**i*t**(2*i+1)/factorial(2*i+1) for i in range(4))))
@lru_cache(None)
def VWF(H):
    b=[Q(0)]+[Q(j+1,N)*CAP[j+1] for j in range(1,M-1)]+[Q(0),Q(1)]
    return tuple(up(sum(comb(j,i)*b[j]*H**(j-i) for j in range(max(1,i),M+1))) for i in range(3))
def verify(outname='reciprocal_strip_certificate.json'):
    co=[None]*N+[Q(1)];co[N-1]=Q(0)
    for j,v in CAP.items():co[j]=v
    var=variance(9,co,U)/2
    assert var<Q(96,100)
    Mtheta=360;Mphase=6;den=512;allh=[]
    for fase in range(Mphase):
        a0=Q(fase,2*N*Mphase);a1=Q(fase+1,2*N*Mphase);hh=[]
        for s in range(Mtheta):
            ta=Q(s,Mtheta);tb=Q(s+1,Mtheta)
            X=sinlower(N*ta,N*tb);Y=sinlower(ta+a0,tb+a1)
            def good(h):
                if not X:return False
                H=Q(h,den);V,W,Z=VWF(H);D=U*Y/N-V
                T=M*H**2*(M/X**2-1)+2*Q(N-2,N)*CAP[N-2]/X
                return D>0 and M*M*D*D>T*(W*W+2*Z*D)
            ll,rr=0,den
            while rr-ll>1:
                mid=(rr+ll)//2
                if good(mid):ll=mid
                else:rr=mid
            assert ll==0 or good(ll)
            hh.append(ll)
        area=2*LOWPI*sum(hh)/(Mtheta*den)
        assert area>Q(94,100),(fase,float(area))
        allh.append(hh)
    print('nonic variance/2',float(var),'six uniform angular area lower bounds',
              [round(float(2*LOWPI*sum(h)/(Mtheta*den)),6) for h in allh])
    (G/outname).write_text(json.dumps({'theta':Mtheta,'phases':Mphase,'den':den,'heights':allh,'variance_display':float(var)},separators=(',',':'))+'\n')
if __name__=='__main__':assert __debug__;verify()

