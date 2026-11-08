"""Independent fixed-grid rational verification. No imports from the first verifier."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import comb,factorial,isqrt
from functools import lru_cache
from pathlib import Path
import argparse,json

UNIT=1<<400
PL=F(3141592653,10**9)
PU=F(3141592654,10**9)
L=F(4,5);C=F(913,1000)
def upper(a):
    a=F(a)
    return F(-((-a.numerator*UNIT)//a.denominator),UNIT)
def sqrt_hi(a):
    a=F(a);m=isqrt(a.numerator*UNIT**2//a.denominator)
    return F(m+1,UNIT)
def small_atan(x,upper_end):
    # Alternating endpoints for nonnegative x<=1/2.
    assert 0<=x<=F(1,2)
    return sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(13 if upper_end else 14))
def angle_hi(x):
    if x<=F(1,2):return small_atan(x,True)
    if x>2:return PU/2-small_atan(1/x,False)
    y=(x-1)/(x+1)
    return PU/4+(small_atan(y,True) if y>=0 else -small_atan(-y,False))
@lru_cache(None)
def angle_loss(H,s):
    # Concavity of arctan on each of the 64 subintervals.
    return upper(sum(upper(angle_hi(H*F(2*j+1,128)/s)) for j in range(64))/(32*PL))
def fixed_row(n):
    if n<=30:return (7,F(1),F(22,100),F(18,100),F(4,100),F(65,64))
    if n==31:return (7,F(104,100),F(17,100),F(16,100),F(4,100),F(65,64))
    if n<=37:return (8,F(96,100),F(15,100),F(16,100),F(15,1000),F(65,64) if n<=36 else F(1))
    k=isqrt(961*n//400)
    k+=400*k*k!=961*n
    return k,F(97,100),F(15,100),F(14,100),F(4,100),F(1)
def check_base():
    assert PL<16*small_atan(F(1,5),False)-4*small_atan(F(1,239),True)
    assert 16*small_atan(F(1,5),True)-4*small_atan(F(1,239),False)<PU
    assert F(433,250)**2<3
    # sum_{p>=0}2^(2p)/(2p+1)! < 11/6, so Cauchy at radius 2 gives |a_l|<=6/2^l.
    terms=[F(4**p,factorial(2*p+1)) for p in range(9)]
    assert sum(terms)+2*terms[-1]*4/(18*19)<F(11,6)
    # Independently recompute the homogeneous coefficient ratio in the energy identity.
    for p in range(16):
        def coeff(j):
            l=p-j
            return (F(comb(2*j,j),(j+1)*4**j)*comb(2*p,2*l)
                  *F(factorial(2*l)*factorial(2*j+2),
                      (4**(p+1))*factorial(l)*factorial(j+1)*factorial(p+1)))
        for j in range(p+1):
            dp=F(comb(2*p+2,p+1),4**(p+1))
            dj=F(comb(2*j+2,j+1),4**(j+1))
            assert coeff(j)*dp==comb(p,j)*coeff(p)*dj
def all_circles(n):
    # Independent geometric partition, using a strictly weaker height estimate.
    lo=F(1,128*n*n);b=F(n-1,2*n);H=sqrt_hi(F(n,3));pieces=0
    while lo<1:
        hi=min(F(1),lo*F(5,4));gain=F(0)
        if hi<=F(1,2):
            gain=6*b*b*max(F(0),b*(1-hi*hi)-H*hi)**2
        need=max(F(0),(F(739,1000)-gain)/L)
        assert need**n<lo,("circle",n,pieces)
        pieces+=1;lo=hi
    return pieces

def run(n,row):
    k,r,s,H,gamma,scale=row;m=n-k
    assert m>=5 and 2<=k<n and min(r,s,H,gamma)>0 and gamma<1 and scale>=0
    y=(PU/(2*k))**2;d=y/6-y*y/120+y**3/5040
    ak=2/PL*d/(1-d);bk=PU/(4*k*k)
    e=upper(2*PU*m*(ak+4*bk*F(95,117)**n/(2*n-1)))
    assert e<2*C
    R=upper(L/(2*C-e))
    S=upper(R*F(1000,433));D=upper((S+R)/(m-1))
    count=m-k-1-(m*angle_loss(H,s)+(2*S+L)/(PL*H))/(1-gamma)
    assert count>0,("real zeros",n,float(count))
    a=scale*F(n-1,m);c=a*R/r;d=a*D/r
    assert d<c<1
    # Alternative exponential lower bound for the product, including the large radius exactly.
    loss=(S-R)*D*(a/r)**2/(2*(1-d));v=c-loss
    assert v>=0
    product=(1-c)*sum(v**i/factorial(i) for i in range(7))
    # Newton recurrence via two running geometric convolutions.
    g=[F(1),F(0)];head=[F(1),F(0)];excess=crit=xx=yy=F(0)
    for q in range(2,m+1):
        xx=upper(R*xx+R*R*g[q-2]);yy=upper(D*yy+(S-R)*D*g[q-2])
        bound=upper(min((xx+yy)/q,F(comb(m,q))*(S/m)**q))
        g.append(bound)
        high=upper(bound*F(comb(n,q),comb(m,q)));head.append(high)
        crit+=F(q,n)*high*r**(n-q)
        excess+=bound*abs(F(comb(n-1,q),comb(m,q))-a**q)/r**q
    diff=product-excess
    assert diff>0,("outer derivative",n,float(diff))
    # Use only a single guaranteed original zero, even if the primary check gets more.
    v=sqrt_hi(s*s+H*H)
    gn=F(0)
    for gj in g:gn=upper(gn*v+gj)
    second=R*R+(S-R)*D
    gn=min(gn,upper(sqrt_hi(s*s+H*H+second/m)**m),upper((v+S/m)**m))
    B=F(0);base=r+H
    for j in range(k):
        q=j;Rolle=F(factorial(n),factorial(m)*factorial(j))
        while q<k:
            h=min(1+q,k-q)
            Rolle*= (s if q==j else 2*s)**h/factorial(h)
            q+=h
        # At most six on |z|=2 in Cauchy's inequality for z/sin(z).
        mult=base**j/(gamma*H)
        for p in range(2,j+1,2):
            mult+=F(6,2**p)*F(factorial(j),factorial(j-p))*H**(p-1)*base**(j-p)
        tail=sum(head[u]*comb(n-u,j)*H**(n-u-j) for u in range(m+1))
        B+=2*(Rolle*gn+tail)*mult
    ratio=B/(n*r**(n-1)*diff);size=crit+B*r*F(n+1,n)
    assert ratio<1 and size<1,("interpolation",n,float(ratio),float(size))
    if 28<=n<=255:
        assert count>F(3,100) and diff>F(19,1000) and ratio<F(9,100) and size<F(71,100)
    return n,float(count),float(diff),float(ratio),float(size),all_circles(n)

if __name__=="__main__":
    arg=argparse.ArgumentParser();arg.add_argument("--low",default=28,type=int)
    arg.add_argument("--high",default=255,type=int);arg.add_argument("--file",default="independent_interpolation_results.json")
    opt=arg.parse_args();assert __debug__;check_base()
    vals=[run(n,fixed_row(n)) for n in range(opt.low,opt.high+1)]
    for i,fn in ((1,min),(2,min),(3,max),(4,max)):
        v=fn(vals,key=lambda a:a[i]);print(i,v[0],v[i],flush=True)
    Path(__file__).with_name(opt.file).write_text(json.dumps(vals,indent=1)+"\n")
    print("Complete consecutive range",opt.low,opt.high,"independently verified.")

