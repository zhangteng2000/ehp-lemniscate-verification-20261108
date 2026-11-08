"""Exact endpoint arithmetic for the all-polynomial line bound in degrees n>=256."""
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
import json

def ceilroot(r,bits=70):
    r=F(r);S=1<<bits
    a=isqrt(r.numerator*S*S//r.denominator)
    return F(a+(a*a*r.denominator!=r.numerator*S*S),S)
def expupper(t):
    t=F(t);assert t>=0
    k=max(40,2*((t.numerator+t.denominator-1)//t.denominator))
    a=s=F(1)
    for i in range(1,k+1):a*=t/i;s+=a
    return s+a*t/(k+1)/(1-t/(k+2))
def expminus(t):
    # An upper bound on exp(-t) from a positive exponential series.
    t=F(t);a=s=F(1)
    for i in range(1,61):a*=t/i;s+=a
    return 1/s
def powupper(base,exponent,bits=55):
    base,exponent=F(base),F(exponent)
    if exponent<0:return powupper(1/base,-exponent,bits)
    p,d=exponent.numerator,exponent.denominator;S=1<<bits
    num=base.numerator**p*S**d;den=base.denominator**p
    lo,hi=0,S
    while hi**d*den<num:hi*=2
    while hi-lo>1:
        a=(lo+hi)//2
        if a**d*den<num:lo=a
        else:hi=a
    return F(hi,S)
def atan_lower_ratio(t,h):
    # For t>=h use atan(t/h)=pi/2-atan(h/t), with the alternating
    # upper polynomial for atan on [0,1].
    z=F(h,t);assert 0<z<=1
    A=sum((-1)**i*z**(2*i+1)/(2*i+1) for i in range(11))
    return 1-2*A/F(333,106)

def verify():
    N=256;u=F(33,256);m=(1-u)*N
    PI=F(355,113);pi0=F(333,106);L=F(4,5);r=F(19,20)
    E=2*PI/4*((PI/12)/(1-PI**2/(24*4*N))+PI/(2*N-1))
    t=(m-2)/2;C=2*t/(1+t+ceilroot(t))
    R=L/(C-E);B=2*(L+E*R)/F(913,1000)
    assert E<F(423,1000) and R<F(577,1000) and B<F(229,100)
    R=F(577,1000);B=F(229,100)
    W=B/(r*(1-u)-R);assert W<F(46,5);W=F(46,5)
    # Split the inverse beta variable at |Z|=1/4 and at lambda=1/2.
    var=u/((1-F(1,N))*(1-u-F(1,N))*N)
    e=F(1,4);assert u<e
    near=W*W/2*expupper(W*e)*var
    prob=expminus(N*e*e*(1-u)**2/(2*(1+e)*(e+u)))
    mid=(expupper(W)+1+W)*prob
    rho=(1+F(1,8))*powupper(8,u)*powupper(F(1,2)+W/m,1-u)
    assert rho<F(9,10) and rho*(N+1)<N
    far=2*N*rho**N
    assert near+mid+far<F(1,2)
    # Upper bound for the high part on every disk of radius rho0..r.
    H=F(1,16);z=F(1,6);outer=F(5,12);rho0=H+outer
    base=2*3*B/(rho0*m)
    tailbase=2*powupper(base,(1-u)/2)
    assert tailbase**N/(1-base)<1 and base<1
    M=expupper((1+2*u/(1-u))*B/rho0)+1
    assert M<500;M=F(500)
    # Oscillation and Blaschke zeros. Positive offsets and a full turn.
    good=F(1,16);A=atan_lower_ratio(z,H)
    frac=(1-u)*A-u-(B/z+3+(2*B+L)/(pi0*H*(1-good)))/N
    assert frac>F(9,25);frac=F(9,25)
    pseud=2*z*outer/(z*z+outer*outer)
    blast=powupper(pseud,frac)*powupper(outer/z,u)
    assert blast<1 and blast**N<F(1,20)
    # Then every good complexified slice of Q is bounded by 2 M rho0^n.
    # Choose two good heights at separation between good*z/(64k)
    # and z/(16k).  The inverse finite difference costs <6/eta.
    spread=(r+H)/(z/2)
    prefactor=2*M*6*64/(good*z)
    ratio=powupper(spread,u)*rho0/r
    assert ratio<1
    assert 4*expupper(W)*prefactor*u*r*ratio**N<1
    # This gives ||Q'|| < n r^(n-1) exp(-W)/4, sufficient for Rouche.
    # A separate upper, needed for critical VALUES, controls the integral
    # of Q' on a radius without relying on the lower top-part bound.
    valratio=powupper(spread,u)*rho0
    assert valratio*(N+1)<N
    value=M*r**N+r*prefactor*u*N*valratio**N
    assert value<F(1,4)
    out=dict(N=N,k="ceil(2 sqrt(n))",curvature_source="all_degree_diagonal_curvature.txt",
        derivative_error_lt="423/1000 times R",root_radius_lt="577/1000",root_sum_lt="229/100",
        W_lt="46/5",beta_remainder_lt="1/2",beta_upper_display=float(near+mid+far),
        beta_parts_display=list(map(float,(near,mid,far))),
        oscillation_fraction_gt="9/25",bl_factor_display=float(blast**N),
        outer_Rouche_ratio_display=float(4*expupper(W)*prefactor*u*r*ratio**N),
        critical_value_radius_lt="1/4",value_bound_display=float(value))
    print(out,flush=True)
    return out
if __name__=="__main__":
    p=Path(__file__).with_name("effective_line_256_certificate.json")
    p.write_text(json.dumps(verify(),indent=2)+"\n")

