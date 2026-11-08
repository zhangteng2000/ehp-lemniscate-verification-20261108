"""Exact strip certificate from the geometric mean of the squared reciprocal roots."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
import json
from verify_reciprocal_strips import sinlower
from verify_interpolation_degrees import LOWPI
from verify_peeled_line import variance
G=Path(__file__).resolve().parent
octcap=[0,0,197,62,319,512,774,0,1000]
def strip_check(n,d,cap,umin,nt,ng,nh):
    m=n-d;top=comb(n,d)
    B=[Q(comb(j+d,d),top)*cap[j+d] for j in range(m+1)]
    @lru_cache(None)
    def coeff(h):
        u=Q(h,nh)
        return tuple(sum(comb(j,l)*B[j]*u**(j-l) for j in range(max(1,l),m+1))
                     for l in range(3))
    # The angular thresholds depend only on |sin(n theta)| and |sin(d theta+gamma)|.
    from math import gcd
    period=Q(gcd(n,d),n)
    records=[]
    for g in range(ng):
        ar=[]
        for i in range(nt):
            a,b=Q(i,nt),Q(i+1,nt);X=sinlower(n*a,n*b)
            Y=sinlower(d*a+g*period/(2*ng),d*b+(g+1)*period/(2*ng))
            def passes(h):
                v,w,z=coeff(h);D=umin*Y/top-v
                return D>0 and m**m*X**2*D**(2*m-2)>(w*w+2*z*D)**m
            l,r=0,nh
            while r-l>1:
                s=(l+r)//2
                if passes(s):l=s
                else:r=s
            assert l==0 or passes(l)
            ar.append(l)
        records.append(ar)
    return records
def main():
    n=8;L=Q(39,50);u=Q(13,25);cap=[Q(v,1000) for v in octcap]
    f=json.loads((G/'octic_energy.json').read_text())['8']
    assert all(Q(f['coefficients'][j])<=cap[j] for j in range(2,n+1))
    v=variance(n,cap,u)/2;assert v<Q(98,100)
    nt,ng,nh=512,8,640
    records=strip_check(n,1,cap,u,nt,ng,nh)
    areas=[2*LOWPI*sum(r)/(nt*nh) for r in records]
    assert min(areas)>L
    out={'n':n,'K':str(L),'u':str(u),'variance_ratio':float(v),'theta':nt,'phases':ng,
         'height_denominator':nh,'strips':records,'area_lower_display':list(map(float,areas))}
    (G/'octic_root_product_strips.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(n,'variance/2',float(v),'certified minimum strip area',float(min(areas)))
if __name__=='__main__':assert __debug__;main()

