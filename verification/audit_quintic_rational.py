"""Independent rational enclosures of trigonometry, special functions, and all quintic analytic constants."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import comb,isqrt
from pathlib import Path
from functools import lru_cache
import json

ROOT=Path(__file__).resolve().parent
P=1<<190

class I:
    __slots__=("l","h")
    def __init__(self,l,h=None):self.l=F(l);self.h=self.l if h is None else F(h);assert self.h>=self.l
    def __add__(self,b):
        b=iv(b);return I(self.l+b.l,self.h+b.h)
    __radd__=__add__
    def __neg__(self):return I(-self.h,-self.l)
    def __sub__(self,b):return self+-iv(b)
    def __rsub__(self,b):return iv(b)+-self
    def __mul__(self,b):
        b=iv(b);v=(self.l*b.l,self.l*b.h,self.h*b.l,self.h*b.h);return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,b):
        b=iv(b);assert b.l*b.h>0;return self*I(1/b.h,1/b.l)
    def __rtruediv__(self,b):return iv(b)/self
    def __pow__(self,k):
        k=F(k)
        if k<0:return I(1)/(self**(-k))
        if k.denominator==1:
            j=k.numerator
            if not j:return I(1)
            v=(self.l**j,self.h**j)
            return I(0 if j%2==0 and self.l<=0<=self.h else min(v),max(v))
        assert self.l>=0
        e=k.denominator
        def nth(t):
            a,b=0,P
            n=t.numerator*P**e;d=t.denominator
            while b**e*d<=n:b*=2
            while b-a>1:
                s=(a+b)//2
                if s**e*d<=n:a=s
                else:b=s
            return a
        l=nth(self.l);h=nth(self.h)
        return I(F(l,P),F(h+1,P))**k.numerator
    def absupper(self):return max(abs(self.l),abs(self.h))
def iv(x):return x if isinstance(x,I) else I(x)

def atan1(n):
    s=sum(F((-1)**k,(2*k+1)*n**(2*k+1)) for k in range(100))
    return I(s,s+F(1,201*n**201))
PI=16*atan1(5)-4*atan1(239)
assert PI.l>F(333,106) and PI.h<F(355,113)


def trig():
    N=4096;S=1<<215;CS=1<<46;oldscale=1<<24
    # sin,cos of pi/2^j from exact half-angle square roots. Never use pi numerically.
    co=(0,0)
    for j in range(2,13):
        old=co
        co=(isqrt((S+old[0])*S//2),1+isqrt(((S+old[1])*S+1)//2))
    si=(isqrt((S-old[1])*S//2),1+isqrt(((S-old[0])*S+1)//2))
    def mul(a,b):
        q=[i*j for i in a for j in b];return min(q)//S,-((-max(q))//S)
    def add(a,b):return a[0]+b[0],a[1]+b[1]
    def neg(a):return -a[1],-a[0]
    x=(S,S);y=(0,0);quarter=[y]
    for j in range(1,N//2):
        x,y=add(mul(x,co),neg(mul(y,si))),add(mul(x,si),mul(y,co))
        quarter.append(y)
    quarter.append((S,S))
    prior=(ROOT/"quintic_sine_dyadic.txt").read_text().splitlines();assert prior[0]=="4096 24"
    high=[]
    for j in range(2*N):
        k=j%N;r=quarter[min(k,N-k)]
        if j>=N:r=neg(r)
        ol,oh=map(int,prior[1+j].split())
        assert ol*S<=r[0]*oldscale<=r[1]*oldscale<=oh*S
        l=r[0]*CS//S;h=-((-r[1]*CS)//S)
        assert h-l<=1
        high.append(str(l)+" "+str(h))
    (ROOT/"quintic_sine_independent.txt").write_text("4096 46\n"+"\n".join(high)+"\n",encoding="ascii")


@lru_cache(None)
def taylor(a,k):
    v=F(1)
    for j in range(k):v*=F(a+j,j+1)
    return v


def J(t):
    n=100;t=F(t)
    def series(a,z):
        v=sum(taylor(a,k)*F(1,2**k)/(2*k+z) for k in range(n))
        return I(v,v+F(2,2**n)/(2*n+z))
    return I(2)**(-t/2)*series(F(1,2),t)+I(2)**F(-1,2)*series(1-t/2,1)


def vin(m,j,l,b):
    a=[F(0)]*(m-j-l+1)
    for k in range(len(a)):a[k]=taylor(j+l-b,k)*taylor(b+j,len(a)-1-k)
    q=[F(0)]*m
    for i,z in enumerate(a):
        for r in range(j):q[l+i+r]+=z*comb(j-1,r)*(-1)**r
    return q


def check():
    x=F(1,5);b=(1-x)/2
    C=I(2)**(1+x)*J(x)
    moment=(I(2)**(2-x))*J(2-x)/PI
    U=2*PI*(1+moment);H=moment**F(1,2)
    assert F(130068,10000)<C.l<C.h<F(130069,10000)
    assert F(137168,10000)<U.l<U.h<F(137169,10000)
    R=I(11);tau=I(F(1,660));z=tau*R
    s=R/(1+z*z)**F(1,2);eps=z/(1-z*z)**F(1,2)
    M=U*(R/2)**x*2*J(1-x)/PI
    core=M*(eps+eps*eps)+(1+x)/(2+x)*tau**2*R**(2+x)
    straight=10*R**(-2+x)/(9*(1-R**-4))**F(1,2)
    v=s**-2;Y=(9/(1-v)+5*v/(1-v)**2)**F(1,2)
    W=max((s**(-2*(j-1)*x)*(j+(5-j)*v)/(1-v)**2).h for j in range(1,6))
    curved=4*J(x)*tau**(1-x)+2*Y*s**(-2+x)+W*s**(-3*x)/(3*x)
    assert (core+straight+curved).h<F(105,100)
    d0=F(4*660,660*660+1); ds=[d0]+[F(1,i) for i in (140,120,96,72,56,44,32,24,16,10,6,4,3)]
    theory=[25,25,25,25,22,16,11,8,5,3,2,1,1]
    # Compute the coefficient polynomials directly over Q and their Toeplitz
    # forms before doing any numerical bound for the gamma quotient.
    sums=[]; normal=[]; K={}
    for m in range(1,26):
        P=[];val=F(0)
        for k in range(m):
            val+=taylor(-b,k)*taylor(b,m-k);P.append(val)
        assert val+taylor(-b,m)==0
        row=[sum(P[i]*taylor(-b,abs(i-k))/taylor(1+b,abs(i-k)) for i in range(m)) for k in range(m)]
        def inner(q):return moment*sum(u*v for u,v in zip(row,q))-P[0]*q[0]
        norm=inner(P);assert norm.l>=0
        normal.append(norm)
        band={}
        for j in range(1,m+1):
            for l in range(m-j+1):
                cost=inner(vin(m,j,l,b)).absupper()
                key=(j,l)
                if key not in K:
                    K[key]=(F(1,5)+1/(4*I(2)**F(1,2)) if key==(1,0) else
                            F((j+l)**(j+l),j**j*(l**l if l else 1))*H)
                band[j+l]=band.get(j+l,F(0))+cost*iv(K[key]).h
        sums.append(band)
    margins=[]
    for i,(d,e,m) in enumerate(zip(ds,ds[1:],theory)):
        h=e/2;ptr=I(0)
        for k in range(m):
            ptr+=(1-h*h)**(k+1)*normal[k]-2*sum(val*h**j for j,val in sums[k].items())
        delta=F(17,16)*I(d/2)**x+2*PI*ptr-(U-C)
        assert delta.l>F(3,1000)
        margins.append(float(delta.l))
    assert (F(17,16)*I(F(1,6))**x-(U-C)).l>F(1,50)
    print("Strict analytic margins from rational enclosures:",min(margins),float((core+straight+curved).h))
    return {"min_chord_margin_gt":"3/1000","small_curvature_loss_lt":"105/100",
            "interval_count":len(margins),"gamma_evaluation":"rational split beta series"}

if __name__=="__main__":
    trig();out=check()
    (ROOT/"quintic_rational_audit.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print("Passed. The dyadic sine table, beta constants and complete chord comparisons were checked.")
