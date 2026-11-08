"""Independent finite arithmetic audit for the threshold proof. Python standard library."""
from fractions import Fraction as Q
from math import comb,isqrt,factorial
import json
from pathlib import Path

def epos(x,M=60):
    x=Q(x);s=sum(x**j/factorial(j) for j in range(M+1))
    return s+x**(M+1)/factorial(M+1)/(1-x/(M+2))
def eneg(x):
    return 1/sum(x**j/factorial(j) for j in range(60))
def rise(x,j):
    p=Q(1)
    for i in range(j):p*= (x+i)/(i+1)
    return p
def lowerroot(q,d):
    S=1<<56;lo,hi=0,2*S
    while hi-lo>1:
        mid=(hi+lo)//2
        if mid**d*q.denominator<=q.numerator*S**d:lo=mid
        else:hi=mid
    return Q(lo,S)
def main():
    N=256;u=Q(33,256);m=Q(223);pc=Q(355,113)
    # Derive both rational pi bounds independently from the Machin identity.
    def atan(x,M):
        p=sum((-1)**j*x**(2*j+1)/(2*j+1) for j in range(M))
        return (p,p+x**(2*M+1)/(2*M+1))
    a,b=atan(Q(1,5),32);c,d=atan(Q(1,239),10)
    assert Q(333,106)<16*a-4*d<=16*b-4*c<pc
    # I = int_0^pi sqrt(sin t) dt. Retain only positive split-beta terms.
    term1=sum(rise(Q(1,2),j)*Q(1,2**j)/(2*j+Q(3,2)) for j in range(18))
    term2=sum(rise(Q(1,4),j)*Q(1,2**j)/(2*j+1) for j in range(18))
    lowI=2*(lowerroot(Q(1,8),4)*term1+lowerroot(Q(1,2),2)*term2)
    assert lowI**2/(2*pc)>Q(913,1000)
    log2lo=2*sum(Q(1,3)**(2*k+1)/(2*k+1) for k in range(7))
    assert 2*pc-8*log2lo<Q(739,1000)
    # Use rounded bounds here, independently of the primary rational routine.
    assert pc/12/(1-pc**2/(24*1024))<Q(263,1000)
    err=pc/2*(Q(263,1000)+pc/511)
    assert err<Q(423,1000)
    assert Q(1052,100)**2>Q(221,2)
    assert 2*Q(221,2)/(1+Q(221,2)+Q(1052,100))>Q(181,100)
    assert Q(8,10)/(Q(181,100)-Q(423,1000))<Q(577,1000)
    assert 2*(Q(4,5)+Q(423,1000)*Q(577,1000))/Q(913,1000)<Q(229,100)
    r=Q(19,20);R=Q(577,1000);B=Q(229,100);W=Q(46,5)
    assert B/(r*(1-u)-R)<W
    var=(Q(33,1)/(255*222))
    near=W**2*epos(W/4)*var/2
    # Use the exact binomial PGF at 1/2 instead of the quadratic Chernoff exponent.
    # log P <= .7(k-2) - (2/5)((n-1)/4 + k-1) <= -(.1-.3u)n.
    mod=(epos(W)+1+W)*eneg(N*(Q(1,10)-Q(3,10)*u))
    # Verify the beta small-lambda exponential directly by taking power 256.
    rho=Q(9,10)
    assert Q(9,8)**256*8**33*(Q(1,2)+W/223)**223<rho**256
    assert 2*N*rho**N+near+mod<Q(1,2)
    assert rho*(N+1)/N<1
    # The angular count uses atan(3/8) by its elementary alternating upper bound.
    z=Q(3,8);at=z-z**3/3+z**5/5-z**7/7+z**9/9
    frac=(1-u)*(1-2*at/Q(333,106))-u-(6*B+3+(2*B+Q(4,5))/(Q(333,106)*Q(1,16)*Q(15,16)))/N
    assert frac>Q(9,25)
    assert (Q(20,29)**(9*256)*Q(5,2)**(33*25))*20**25<1
    ro=Q(23,48)
    assert epos((1+2*u/(1-u))*B/ro)+1<500
    tail=2*3*B/(ro*m)
    assert tail<1 and 2**256*tail**111/(1-tail)<1
    # For n>=256 the tail base is <1, hence testing its endpoint suffices.
    assert (2**256*tail**111)<1
    spread=(Q(19,20)+Q(1,16))*12
    assert (ro/r)**256*spread**33<Q(3,4)**256
    pref=Q(1000)*6*64*16*6
    assert 4*epos(W)*pref*u*r*Q(3,4)**N<1
    assert 500*r**N+r*pref*u*N*(r*Q(3,4))**N<Q(1,4)
    assert (r*Q(3,4))*(N+1)<N
    assert Q(4,5)**256*Q(15,16)**256>Q(739,1000)**256
    assert Q(15,16)**256 < Q(1,128*256**2)
    # Check a few arbitrary-degree binomial integral identities for beta inverse moments.
    for n,k in [(13,4),(17,7),(28,5)]:
        m=n-k
        for j in range(m+1):
            assert Q(comb(n-1,j),comb(m,j))==Q(factorial(m-j)*factorial(n-1),factorial(m)*factorial(n-1-j))
        lam=Q(n-1,m)
        assert Q(comb(n-1,2),comb(m,2))/lam**2-1==Q(k-1,(m-1)*(n-1))
    print("Independent analytic endpoint checks passed, including c*, root budgets, beta moments, oscillation, Rouche and the circle gap.")
    print("Alternative beta error <",float(near+mod+2*N*rho**N))
    Path(__file__).with_name("effective_line_256_independent.json").write_text(json.dumps({
      "n0":256,"alternate_binomial_bound":str(Q(1,2)),"beta_upper_display":float(near+mod+2*N*rho**N),
      "matching_constant_gt":"913/1000","oscillation_fraction_gt":"9/25","source_geometry":True},indent=2)+"\n")
if __name__=="__main__":main()

