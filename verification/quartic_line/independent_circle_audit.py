"""A second quartic circle comparison, using coarse rational constants.

This file imports no verification helper.  It deliberately uses larger
constants and different geometric-series tails than the primary curvature
check.  All inequalities used for acceptance are exact rational inequalities.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import isqrt,comb
from pathlib import Path
import json

S=1<<96

def root(q,d):
    q=F(q);assert q>=0 and d in (2,4,8)
    z=q.numerator*S**d//q.denominator
    a=z
    for _ in range(d.bit_length()-1):a=isqrt(a)
    assert F(a**d,S**d)<=q<F((a+1)**d,S**d)
    return F(a,S),F(a if F(a**d,S**d)==q else a+1,S)

def positive_series(a,z,N):
    # sum (a)_j 2^-j /(j! (z+j)), with 0<a<=1.
    assert 0<a<=1 and z>0
    v=F(1);s=F(0)
    for j in range(N):
        s+=v/(z+j)
        v*=F(a+j,2*(j+1))
    return s,s+F(2,2**N)/(z+N)

def beta(a,N=40):
    l1,h1=positive_series(F(1,2),a,N)
    l2,h2=positive_series(1-a,F(1,2),N)
    # Only a=1/8 or 3/8 is needed here.
    d=8;ll,hh=root(2,d)
    power=int(8*a);slo,shi=1/hh**power,1/ll**power
    rlo,rhi=root(2,2)
    return slo*l1+l2/rhi,shi*h1+h2/rlo

def verify():
    pi_lo,pi_hi=F(333,106),F(355,113)
    # Machin's formula with alternating rational remainders proves the pi bounds.
    def atan(q,N):
        s=sum(F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N))
        e=F(1,(2*N+1)*q**(2*N+1))
        return (s,s+e) if N%2==0 else (s-e,s)
    a,b=atan(5,44);c,d=atan(239,12)
    assert 16*a-4*d>pi_lo and 16*b-4*c<pi_hi
    assert beta(F(3,8))[1]<F(39,10)
    assert beta(F(1,8))[1]<F(19,2)
    # Independent upper bound for U4 from the positive binomial square series.
    b=F(3,8);q=F(1);moment=F(1);N=160
    for j in range(1,N+1):
        q*=F(b-j+1,j);moment+=q*q
    moment+=N*q*q  # The standard k^-2 tail, proved in the analytic note.
    U=2*pi_hi*(1+F(3,4)*moment)
    assert U<F(2353,200)
    C=root(2,4)[0]*beta(F(1,8))[0]
    assert C>F(1107,100)
    assert U-C<F(139,200)

    x=F(1,4);R=F(5,2);tau=F(1,64);rc=F(5,3)
    H2=F(7,8);H=F(15,16);s=F(249,100)
    assert H*H>H2 and s*s*(1+tau*tau*R*R)<R*R
    assert U/(2*pi_lo)-1<H2
    e=(tau*R)/root(1-tau*tau*R*R,2)[0]
    M=F(2353,200)*root(R/2,4)[1]*F(39,10)/pi_lo
    core=M*(e+e*e)+F(5,9)*tau*tau*R*R*root(R,4)[1]
    ratio=(rc/R)**4
    linear=H*root(R,4)[1]*(rc/R)**2/(F(7,4)*root(1-ratio,2)[0])
    quadratic=2*H2*rc/(F(3,4)*root(R,4)[0]**3*(1-ratio))
    motion=(4*linear+2*quadratic)*tau/(F(2,5)-tau)
    first=2*F(35,512)*F(19,2)*root(tau,4)[1]**3
    def K(v,l):
        assert v>1
        return 1+1/v+(1+pi_hi*l/2)/min(v-1,1)
    y=1/root(s,4)[0]
    top=200
    rows=[]
    for j in range(2,top):
        val=F(0)
        for k in range(2,top):
            if (j-k)%4 or j+k<=5:continue
            val+=K(F(j+k-1,4),F(abs(j-k),4))*y**(j+k)
        rows.append(val)
    # K((j+k-1)/4,|j-k|/4) <= 6+2(j+k).
    # Bound the sum of EVERY omitted matrix entry by the sum over j+k>=top+2.
    m=top+2
    geom=1/(1-y)
    hsum=y**m*(m*geom+y*geom**2)
    hsq=y**m*(m*m*geom+2*m*y*geom**2+y*(1+y)*geom**3)
    row=max(rows)+6*hsum+2*hsq
    Q=24;z=1/s
    square=sum((K(F(q)-x,q)*z**q)**2 for q in range(2,Q))
    tail=z**Q*((3+2*Q)/(1-z)+2*z/(1-z)**2)
    A2=root(square,2)[1]+tail
    higher=2*tau*R*root(R,4)[1]*(H2*row+2*H*A2)
    pieces=dict(core=core,motion=motion,first=first,higher=higher)
    total=sum(pieces.values())
    limits={'core':F(631,1000),'motion':F(253,1000),'first':F(58,1000),'higher':F(433,1000)}
    assert all(pieces[k]<limits[k] for k in pieces)
    assert total<F(7,5)<F(29,20)
    delta=F(128,4097);g=F(12157,32776)
    assert F(21,50)**4<delta
    margin=2*pi_lo*F(3,8)**2*g*g+F(29,20)*F(21,50)-F(139,200)
    assert margin>F(1,32)
    assert F(29,80)>2*pi_hi*F(3,8)**3
    assert (F(29,20)-F(7,5))/F(6,5)>F(1,32)
    assert (F(29,24)-F(139,200))/F(6,5)>F(1,32)
    def upper(q):return F((q*(1<<40)).__ceil__(),1<<40)
    out={'status':'all independent rational checks passed',
         'curvature_component_upper':{k:str(upper(v)) for k,v in pieces.items()},
         'curvature_total_upper':str(upper(total)),
         'curvature_bound':'7/5','quartic_gap_constant':'1/32',
         'U4_upper':'2353/200','C4_lower':'1107/100',
         'matrix_cutoff':top,'linear_cutoff':Q}
    Path(__file__).with_name('independent_circle_audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Independent curvature upper:',float(total),' < 7/5; all scalar checks passed.')
    return out

if __name__=='__main__':assert __debug__;verify()
