"""Exact verification of the interpolation criterion, including the circle transfer."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial,isqrt
from pathlib import Path
import argparse
import json

HERE=Path(__file__).resolve().parent
LOWPI=Q(3141592653589793238,10**18)
HIGHPI=LOWPI+Q(1,10**18)
L=Q(4,5)
C=Q(913,1000)

def up(x,bits=66):
    x=Q(x)
    assert x>=0
    if not x:return x
    shift=bits+x.denominator.bit_length()-x.numerator.bit_length()
    unit=Q(1,1<<shift) if shift>=0 else Q(1<<(-shift))
    return -(-(x/unit).numerator//(x/unit).denominator)*unit

def sqrt_up(x):
    x=Q(x)
    if not x:return x
    n=isqrt((x.numerator<<160)//x.denominator)
    return Q(n+(n*n*x.denominator!=(x.numerator<<160)),1<<80)

def atansum(x,K):
    return sum((-1)**i*x**(2*i+1)/(2*i+1) for i in range(K+1))

def pi_check():
    lo=16*atansum(Q(1,5),17)-4*atansum(Q(1,239),6)
    hi=16*atansum(Q(1,5),18)-4*atansum(Q(1,239),7)
    assert LOWPI<lo<hi<HIGHPI
    assert Q(433,250)**2<3

def atan_up(x):
    if x>2:return HIGHPI/2-atansum(1/x,17)
    if x<=Q(1,2):return atansum(x,18)
    u=(x-1)/(x+1)
    return HIGHPI/4+(atansum(u,18) if u>=0 else -atansum(-u,17))

@lru_cache(None)
def chi_up(x):
    v=x*x/(2+x*x)
    half_log=sum(v**(2*i+1)/(2*i+1) for i in range(50))
    return up(2/LOWPI*(atan_up(x)-half_log/x))

def constants(n,k):
    m=n-k
    v=LOWPI/(2*k)
    sinlo=sum((-1)**i*v**(2*i+1)/factorial(2*i+1) for i in range(6))
    ak=1/(k*sinlo)-2/HIGHPI
    bk=HIGHPI/(4*k*k)
    err=up(2*HIGHPI*m*(ak+4*bk*Q(95,117)**n/(2*n-1)))
    assert err<2*C and m>=5
    R=up(L/(2*C-err));S=up(Q(1000,433)*R);D=up((S+R)/(m-1))
    assert D<R<S
    # The explicit fractional moment bounds remain valid if rounded upward.
    moments=[0]+[up(R**j+(S-R)*D**(j-1)) for j in range(1,m+1)]
    g=[Q(1),Q(0)]
    powers=[Q(1)]
    for j in range(1,m+1):
        powers.append(up(powers[-1]*S/m))
        if j>=2:
            newton=sum(moments[h]*g[j-h] for h in range(2,j+1))/j
            g.append(up(min(comb(m,j)*powers[j],newton)))
    head=[up(g[j]*Q(comb(n,j),comb(m,j))) for j in range(m+1)]
    return R,S,D,moments,g,head

@lru_cache(None)
def inverse_sinc(k):
    c=[Q(1)]
    for j in range(1,k):
        c.append(-sum(Q((-1)**h,factorial(2*h+1))*c[j-2*h] for h in range(1,j//2+1)))
    return c

def circle(n):
    start=Q(1,128*n*n);b=Q(n-1,2*n);K=sqrt_up(n/LOWPI)
    # Curvature at tau_0: 1/2 + 8/128^(5/6) < 4/5.
    assert (Q(8)/(L-Q(1,2)))**6<128**5
    stack=[(start,Q(1,2)),(Q(1,2),Q(1))];leaves=0
    while stack:
        lo,hi=stack.pop()
        extra=Q(0)
        if hi<=Q(1,2):
            extra=2*LOWPI*b*b*max(Q(0),b*(1-hi*hi)-K*hi)**2
        q=(Q(739,1000)-extra)/L
        if q<=0 or lo>q**n:
            leaves+=1
        else:
            assert hi-lo>Q(1,2**28),("circle",n,lo,hi)
            h=(lo+hi)/2;stack.extend([(lo,h),(h,hi)])
    return leaves

def verify(n,row):
    k=int(row[0]);r,s,H,gamma,scale=map(Q,row[1:])
    m=n-k
    assert 2<=k<=n-5 and 0<gamma<1 and min(s,r,H)>0 and scale>=0
    R,S,D,P,g,T=constants(n,k)
    count=m-k-1-(m*chi_up(H/s)+(2*S+L)/(LOWPI*H))/(1-gamma)
    assert count>0,("count",n,float(count))
    M=min(k,-((-count.numerator)//count.denominator))
    # The simple certificate can always choose M=1. More zeros only improve it.
    v=sqrt_up(s*s+H*H)
    gp=min(up(sum(g[i]*v**(m-i) for i in range(m+1))),
           up((v+S/m)**m),up(sqrt_up(s*s+H*H+P[2]/m)**m))
    pcoef=[]
    for j in range(k):
        i=j;blocks=[]
        while i<k:
            size=min(i+M,k-i);blocks.append(size);i+=size
        den=1
        for size in blocks:den*=factorial(size)
        factor=Q(comb(n,k)*comb(k,j)*factorial(k-j)*2**(k-j-blocks[0]),den)
        real=gp*factor*s**(k-j)
        high=sum(T[q]*comb(n-q,j)*H**(n-q-j) for q in range(m+1))
        pcoef.append(up(real+high))
    rad=r+H;alpha=inverse_sinc(k)
    V=up(2*sum(pcoef[j]*(rad**j/(gamma*H)+sum(abs(alpha[q])*comb(j,q)*factorial(q)*
             H**(q-1)*rad**(j-q) for q in range(2,j+1,2))) for j in range(k)))
    a=scale*Q(n-1,m);c1=a*R/r;c2=a*D/r
    assert max(c1,c2)<1,("outer",n)
    x=c1-(S-R)*D*(a/r)**2/(2*(1-c2))
    prod=(1-c1)*sum(x**j/factorial(j) for j in range(8))
    error=sum(g[j]*abs(Q(comb(n-1,j),comb(m,j))-a**j)/r**j for j in range(2,m+1))
    lower=prod-error
    ratio=V/(n*r**(n-1)*lower) if lower>0 else Q(10)
    val=sum(Q(j,n)*T[j]*r**(n-j) for j in range(2,m+1))+r*(1+Q(1,n))*V
    assert lower>0 and ratio<1 and val<1,("Rouche/critical values",n,float(lower),float(ratio),float(val))
    if 28<=n<=255:
        assert count>Q(3,100) and lower>Q(19,1000) and ratio<Q(9,100) and val<Q(71,100)
    return dict(n=n,interior_zeros=M,raw_count_lower_display=float(count),
        normalized_T_derivative_lower_display=float(lower),Rouche_ratio_display=float(ratio),
        critical_radius_upper_display=float(val),circle_intervals=circle(n))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--params",default="interpolation_degree_witnesses.json")
    p.add_argument("--out",default="interpolation_degree_audit.json")
    p.add_argument("--low",type=int,default=28);p.add_argument("--high",type=int,default=255)
    u=p.parse_args();assert __debug__
    pi_check();rows=json.loads((HERE/u.params).read_text());results=[]
    for n in range(u.low,u.high+1):
        results.append(verify(n,rows[str(n)]))
    assert len(results)==u.high-u.low+1
    for name,op in (("raw_count_lower_display",min),("normalized_T_derivative_lower_display",min),
                    ("Rouche_ratio_display",max),("critical_radius_upper_display",max)):
        w=op(results,key=lambda row:row[name]);print(name,w["n"],w[name],flush=True)
    (HERE/u.out).write_text(json.dumps(results,indent=1)+"\n")
    print("All degrees",u.low,"through",u.high,"passed",flush=True)

