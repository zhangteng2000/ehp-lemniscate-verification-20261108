"""Rational scalar certificate for collinear critical values by successive differentiation."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from math import comb,factorial
from pathlib import Path
from functools import lru_cache
import argparse,json
from verify_interpolation_degrees import up,sqrt_up,LOWPI,HIGHPI,pi_check
HERE=Path(__file__).resolve().parent
C=Q(913,1000);NU=Q(1000,433);L=Q(49,60)
def choose(n,j):return comb(n,j) if 0<=j<=n else 0
def vmax(m,j,tail,amin,u):
    return sum(tail[q]*amin**q*(u/(m-1))**(j-q)*
       abs(choose(m-1-q,j-q)-(m-1)*choose(m-1-q,j-q-1)) for q in range(min(j,m-1)+1))
def slopes(cap):
    hull=[]
    for i,v in enumerate(cap):
        hull.append((i,v))
        while len(hull)>=3:
            a,b,c=hull[-3:]
            if (b[1]-a[1])*(c[0]-b[0]) >= (c[1]-b[1])*(b[0]-a[0]):break
            hull.pop(-2)
    lam=[]
    for a,b in zip(hull,hull[1:]):
        lam += [(b[1]-a[1])/(b[0]-a[0])]*(b[0]-a[0])
    assert all(a>=b>=0 for a,b in zip(lam,lam[1:]))
    for i,v in enumerate(cap):assert sum(lam[:i])>=v
    return lam

@lru_cache(None)
def envelope(m):
    if m==1:return (Q(1),Q(0))
    if m==2:return (Q(1),Q(0),Q(1))
    den=24;tail=envelope(m-1);upper=[Q(1),Q(0)]+[Q(0)]*(m-1)
    weights={Q(0),Q(1,4),Q(1,2),Q(1),Q(2),Q(3),Q(4),Q(5),Q(7),Q(9),Q(12),Q(16),
             Q(m-1,4),Q(m-1,2),Q(3*(m-1),4)}
    for i in range(den):
        lo=Q(i,den);hi=Q(i+1,den);cap=[Q(0)]
        for h in range(1,m):
            v=min(h*hi,NU-lo)
            if h<m-1:v=min(v,1+Q(2*h,m-h-1)-lo)
            if m>=5:v=min(v,1+h*(NU+1)/(m-1)-lo)
            cap.append(v)
        lam=slopes(cap)
        powers=[0]+[up(sum(v**j for v in lam)) for j in range(1,m+1)]
        expo=[Q(1),Q(0)]
        for j in range(2,m+1):expo.append(up(sum(powers[q]*expo[j-q] for q in range(2,j+1))/j))
        phase=[Q(1),Q(0)]+[Q(j-1,factorial(j))*hi**j for j in range(2,m+1)]
        inherited=min((1-Q(m,m-1)*lo*t/(1+t))/(1-t/(m-1)) for t in weights if t<m-1)
        assert inherited>0
        for j in range(2,m+1):
            convol=sum(phase[q]*expo[j-q] for q in range(j+1))
            moved=vmax(m,j,tail,inherited,hi)
            size=up(min(convol,moved,comb(m,j)*min(NU/m,hi)**j))
            upper[j]=max(upper[j],size)
    return tuple(upper)

@lru_cache(None)
def Ebound(n,k):
    m=n-k;x=LOWPI/(2*k)
    bot=sum((-1)**i*x**(2*i+1)/factorial(2*i+1) for i in range(6))
    y=HIGHPI/(2*k)
    coslo=sum((-1)**i*y**(2*i)/factorial(2*i) for i in range(6))
    sinhi=sum((-1)**i*y**(2*i+1)/factorial(2*i+1) for i in range(7))
    bk=2/LOWPI-coslo/(k*sinhi)
    return up(2*HIGHPI*m*(1/(k*bot)-2/HIGHPI+4*bk*Q(95,117)**n/(2*n-1)))

def caps(n,wit=None):
    co=[None]*(n-1)+[Q(0),Q(1)]
    stored=[]
    for lap in range(4):
        for k in range(n-2,1,-1):
            m=n-k;E=Ebound(n,k);e=envelope(m)
            prior=[None if co[n-j] is None else co[n-j]*Q(comb(m,j),comb(n,j)) for j in range(m+1)]
            def test(r):
                a=(L/r+E)/(2*C)
                return sum(e[j]*a**j if prior[j] is None else min(prior[j]/r**j,e[j]*a**j)
                           for j in range(2,m+1))<1
            if wit is None:
                base=L/(2*C-E) if E<2*C else None
                hi=Q(1)
                while not test(hi):
                    hi*=2
                    if hi>=64:break
                if test(hi):
                    lo=Q(0)
                    for _ in range(22):
                        mid=(lo+hi)/2
                        if test(mid):hi=mid
                        else:lo=mid
                    r=hi if base is None else min(hi,base)
                else:
                    assert base is not None,(n,k)
                    r=base
                r=Q((r*2**22).__ceil__(),2**22)
            else:
                r=Q(wit[len(stored)])
            assert r>0 and ((E<2*C and r*(2*C-E)>=L) or test(r)),(n,k)
            stored.append(str(r))
            A=up((L+E*r)/(2*C))
            for j in range(2,m+1):
                val=up(e[j]*A**j*Q(comb(n,j),comb(m,j)))
                co[n-j]=val if co[n-j] is None else min(co[n-j],val)
    assert all(x is not None for x in co[2:])
    return co,stored

def moments(d,N):
    m=len(d)-1;p=[Q(0)]
    for q in range(1,N+1):
        p.append(sum(d[i]*(i if i==q else p[q-i]) for i in range(1,min(m,q)+1)))
    return p

def variance(n,co,u):
    a=co.copy();a[1]=u;m=n-1
    ds=[0]+[Q(n-i,n)*a[n-i] for i in range(1,n)]
    p=moments(ds,2*n-2)
    out=sum(Q(n-i,n)*Q(n-j,n)*a[i]*a[j]*(p[i+j]-p[i]*p[j]/m)
            for i in range(1,n-1) for j in range(1,n-1))
    assert out>=0
    return out

def root_below(v,m):
    S=1<<50;lo,hi=0,2*S
    while hi**m*v.denominator<=v.numerator*S**m:hi*=2
    while hi-lo>1:
        q=(hi+lo)//2
        if q**m*v.denominator<=v.numerator*S**m:lo=q
        else:hi=q
    return Q(lo,S)

def gram(n,co,u,rad):
    m=n-1;t=root_below(u/n,m)
    assert 0<t and t**m<=u/n and rad>1
    vals=[0,0]+[up(co[j]*t**(j-n)) for j in range(2,n-1)]
    d=[0,Q(0)]+[up(Q(n-q,n)*vals[n-q]) for q in range(2,n-1)]+[Q(1)]
    assert len(d)==m+1
    assert sum(d[i]/rad**i for i in range(1,n))<1,("outer radius",n,float(rad))
    p=moments(d,m);B=sqrt_up(sum(x*x for x in d[1:-1]))
    need=m*(m-p[2]);assert need>0
    used=sum(Q(n-j,n)*vals[j]*(p[j-1]+p[j+1]+B*(1+(rad+1)**2*rad**(j-1)/4))
             for j in range(2,n-1))
    return used/need

def run(opts):
    pi_check()
    pars=json.loads((HERE/opts.params).read_text())
    w={};rows=[]
    for n in range(opts.low,opts.high+1):
        u,r=map(Q,pars[str(n)]);c,ws=caps(n)
        sm=variance(n,c,u)/2;lg=gram(n,c,u,r)
        assert sm<Q(94,100) and lg<Q(84,100),("line",n,float(sm),float(lg))
        print(n,"variance/2",float(sm),"Gram ratio",float(lg),flush=True)
        w[str(n)]=ws;rows.append([n,float(sm),float(lg)])
    (HERE/opts.prefix).with_suffix(".json").write_text(json.dumps({"radii":w,"checks":rows},separators=(',',':'))+"\n")
if __name__=="__main__":
    pg=argparse.ArgumentParser();pg.add_argument("--low",type=int,default=15);pg.add_argument("--high",type=int,default=27)
    pg.add_argument("--params",default="collinear_gram_parameters.json");pg.add_argument("--prefix",default="peeled_primary")
    o=pg.parse_args();assert __debug__;run(o)

