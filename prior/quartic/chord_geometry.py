"""Certified quartic shape charts over the triangle of three target values."""
from fractions import Fraction as F
from pathlib import Path
import sys
import os
import math

vendor=os.environ.get("EHP_PYTHON_FLINT")
if vendor:sys.path.insert(0,vendor)
else:
    candidate=Path(__file__).resolve().parent.parent.parent/"research_20261105"/"vendor_flint"
    if candidate.exists():sys.path.insert(0,str(candidate))
from flint import arb,acb,ctx
import slice_data as s
ctx.prec=104
OMEGA=acb(-arb(1)/2,arb(3).sqrt()/2)
I=acb(0,1)
ALPHA=arb(5)/2


class Split(Exception):pass
class Outside(Exception):pass


def real(x):
    return arb(x.numerator)/x.denominator


def interval(x,h):
    return arb(real(x),real(h))


def rational(x):
    m,e=x.man_exp()
    return F(int(m))*F(2)**int(e)


def fixed_down(x):
    return int(x.lower().floor().unique_fmpz())


def fixed_up(x):
    return int(x.upper().ceil().unique_fmpz())


class J:
    __slots__=("v","s","t")
    def __init__(self,v,ds=0,dt=0):
        self.v=acb(v);self.s=acb(ds);self.t=acb(dt)
    def __add__(self,b):
        if isinstance(b,J):return J(self.v+b.v,self.s+b.s,self.t+b.t)
        return J(self.v+b,self.s,self.t)
    __radd__=__add__
    def __neg__(self):return J(-self.v,-self.s,-self.t)
    def __sub__(self,b):return self+(-b)
    def __rsub__(self,b):return -self+b
    def __mul__(self,b):
        if isinstance(b,J):
            return J(self.v*b.v,self.s*b.v+self.v*b.s,self.t*b.v+self.v*b.t)
        return J(self.v*b,self.s*b,self.t*b)
    __rmul__=__mul__
    def inv(self):
        q=1/self.v
        if not q.is_finite():raise Split("division")
        return J(q,-self.s*q*q,-self.t*q*q)
    def __truediv__(self,b):return self*b.inv() if isinstance(b,J) else self*(1/arb(b))
    def conj(self):return J(self.v.conjugate(),self.s.conjugate(),self.t.conjugate())
    def re(self):return J(self.v.real,self.s.real,self.t.real)
    def im(self):return J(self.v.imag,self.s.imag,self.t.imag)
    def norm2(self):return (self*self.conj()).re()
    def root(self,n):
        if not self.v.real>0:raise Split("radical branch")
        q=self.v.root(n);k=q/(n*self.v)
        return J(q,k*self.s,k*self.t)


def phase(sv,tv):
    h=(1+tv*I)/(1-tv*I)
    w=h*h*h;r=sv*sv*sv
    return h,w,r,r*w


def algebra(sh,sv,tv):
    h,w,r,v=phase(sv,tv)
    g=sv*h*((1-v)*4).root(3)
    y=(1-g*OMEGA).root(2);z=(1-g*OMEGA.conjugate()).root(2)
    x=(1-v*2)/(y*z)
    if sh==0:c=x+y+z
    elif sh==1:c=x-y-z
    elif sh==2:c=-x+y-z
    else:c=-x-y+z
    d=c*F16
    sin=w.im();m=(d*(sin+(r-w.re())*I))*(1/ALPHA**4)
    mc=m.conj()*c;norm=c.norm2()
    vec=[sin,2*sin*c.re(),-2*sin*c.im(),sin*norm,-m.re(),-m.im(),-mc.re(),mc.im()]
    clip=[J(1),2*c.re(),-2*c.im(),norm,-d.norm2()*(arb(s.T)**2/ALPHA**8)]
    return vec,clip,c,v,w

F16=arb(16)/3


def geometry(sh,ss,tt,hs,ht):
    sv=J(interval(ss,hs),1,0);tv=J(interval(tt,ht),0,1)
    h,w,r,v=phase(sv,tv)
    if v.v.real>arb(1)/2 or w.v.real<0:
        raise Outside
    # |1-v|<=1 is equivalent to rho<=2 cos(theta), including the limit rho=0.
    if r.v.real>2*w.v.real:raise Outside
    b0,clips0,c0,v0,w0=algebra(sh,J(real(ss)),J(real(tt)))
    bv,clips,c,v,w=algebra(sh,sv,tv)
    ds,dt=real(hs),real(ht)
    # For scalar geometric cutoffs use both direct and differential enclosures.
    def tight(f,p):
        width=abs(f.s)*ds+abs(f.t)*dt
        return p.v+acb(arb(0,width.upper()),arb(0,width.upper()))
    def mag(f,p):
        a=abs(f.v);b=abs(tight(f,p))
        return max(F(0),rational(a.lower()),rational(b.lower())),min(rational(a.upper()),rational(b.upper()))
    sv0=J(real(ss));tv0=J(real(tt))
    wval=w.im();wc=w0.im()
    sinlow=max(F(0),rational(wval.v.real.lower()),rational(tight(wval,wc).real.lower()))
    sinup=min(F(1),rational(wval.v.real.upper()),rational(tight(wval,wc).real.upper()))
    dl,du=mag(1-v,1-v0)
    dl=max(dl,F(1,2));du=min(du,F(1))
    if du<dl or sinup<sinlow:raise Outside
    klo=2*sinlow/du;khi=2*sinup/dl
    cl,cu=mag(c,c0)
    # The depressed quartic for c gives 1/6<|c|<4 on the admissible lens.
    cl=max(cl,F(1,6));dslo=cl*16/3
    full=klo>=F(1,80)
    if full:
        target=s.MODEL_LOWER*s.fourthlo(1/khi)
    else:
        q=min(F(1),s.T*klo/2)
        # For performance, round the lower input down before the powers.
        q=F((q*(1<<30)).__floor__(),1<<30);zz=q*q
        val=F(0)
        for a in reversed(s.TWEIGHTS):val=val*zz+a
        target=24*F((val*(1<<36)).__floor__(),1<<36)
        if khi*s.T>=2:target=min(target,s.MODEL_LOWER*s.fourthlo(1/khi))
        target-=F(40,81)*s.sqrtup(1-zz)
    target*=s.fourthlo(dslo)*F(2,5)
    def emit(vec,p):
        C=1<<s.SHAPE_BITS
        norm=max(float(abs(k.v.real).upper()) for k in p)
        if not math.isfinite(norm) or norm==0:raise Split("shape scale")
        _,ex=math.frexp(norm)
        scale=arb(2)**(s.SHAPE_BITS-ex)
        gens=[[],[],[]];remainder=[]
        for f,f0 in zip(vec,p):
            terms=[f0.v.real*scale,f.s.real*ds*scale,f.t.real*dt*scale]
            err=arb(0)
            for i,z in enumerate(terms):
                a=fixed_down(z.mid());gens[i].append(a)
                err+=abs(z-a)
            remainder.append(fixed_up(err))
        if sum(sum(map(abs,row)) for row in gens)+sum(remainder)>1<<32:raise Split("large shape interval")
        return gens,remainder
    gens,rad=emit(bv,b0)
    cc=[];cr=[]
    for f,p in zip(clips,clips0):
        a=p.v.real
        err=abs(f.s.real)*ds+abs(f.t.real)*dt
        fac=1<<12
        l=max(fixed_down(f.v.real*fac),fixed_down((a-err)*fac))
        h=min(fixed_up(f.v.real*fac),fixed_up((a+err)*fac))
        cval=(l+h)//2
        cc.append(cval);cr.append(max(cval-l,h-cval))
    if sum(abs(a)+b for a,b in zip(cc,cr))>1<<32:raise Split("large clip interval")
    # Heuristic anisotropic split scores. Only the enclosure above is used to prove length.
    src=(5000,1600,1600,500,60,60,20,20) if not full else (80,40,40,16,10,10,5,5)
    scores=[sum(abs(a)*t for a,t in zip(row,src)) for row in gens[1:]]
    qerr=F(target.numerator*(1<<38)//target.denominator,1<<38)
    return full,qerr,gens,rad,cc,cr,scores,klo,khi

