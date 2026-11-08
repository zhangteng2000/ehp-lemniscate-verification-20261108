"""Outward dyadic real/complex intervals. Integer arithmetic only."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import isqrt
BITS=60
S=1<<BITS

class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None,raw=False):
        if raw:self.lo,self.hi=int(lo),int(hi)
        else:
            lo=F(lo);hi=lo if hi is None else F(hi)
            self.lo=(lo*S).__floor__();self.hi=(hi*S).__ceil__()
        assert self.lo<=self.hi
    def __add__(self,b):
        b=iv(b);return I(self.lo+b.lo,self.hi+b.hi,raw=True)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo,raw=True)
    def __sub__(self,b):return self+-iv(b)
    def __rsub__(self,b):return iv(b)+-self
    def __mul__(self,b):
        b=iv(b);v=[x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return I(min(v)//S,-((-max(v))//S),raw=True)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=iv(b);assert b.lo*b.hi>0
        if b.hi<0:return (-self)/(-b)
        return self*I(S*S//b.hi,-((-S*S)//b.lo),raw=True)
    def __rtruediv__(self,b):return iv(b)/self
    def sq(self):
        vals=(self.lo*self.lo,self.hi*self.hi)
        return I(0 if self.lo<=0<=self.hi else min(vals)//S,-((-max(vals))//S),raw=True)
    def sqrt(self):
        assert self.lo>=0
        a=isqrt(self.lo*S);b=isqrt(self.hi*S)
        if b*b<self.hi*S:b+=1
        return I(a,b,raw=True)
    def bounds(self,bits):
        assert bits<=BITS
        d=1<<(BITS-bits);return self.lo//d,-((-self.hi)//d)
    def lower(self):return F(self.lo,S)
    def upper(self):return F(self.hi,S)

def iv(a):return a if isinstance(a,I) else I(a)

class C:
    __slots__=('re','im')
    def __init__(self,re=0,im=0):self.re,self.im=iv(re),iv(im)
    def __add__(self,b):
        b=cv(b);return C(self.re+b.re,self.im+b.im)
    __radd__=__add__
    def __neg__(self):return C(-self.re,-self.im)
    def __sub__(self,b):return self+-cv(b)
    def __rsub__(self,b):return cv(b)+-self
    def __mul__(self,b):
        b=cv(b);return C(self.re*b.re-self.im*b.im,self.re*b.im+self.im*b.re)
    __rmul__=__mul__
    def __truediv__(self,b):
        if isinstance(b,C):return self*b.conj()/b.abs().sq()
        return C(self.re/b,self.im/b)
    def conj(self):return C(self.re,-self.im)
    def abs(self):return (self.re.sq()+self.im.sq()).sqrt()

def cv(a):return a if isinstance(a,C) else C(a)

def shape(lo,hi,sheet):
    """All c in one of the four closed real-critical-value charts."""
    t=I(lo,hi);g=1-t.sq();rr=(1+g+g.sq()).sqrt()
    p=((rr+1+g/2)/2).sqrt();q=I(3).sqrt()*g/(2*p)
    if sheet==0:c=C(t+2*p)
    elif sheet==1:c=C(t-2*p)
    elif sheet==2:c=C(-t,q)
    elif sheet==3:c=C(-t,-q)
    else:raise ValueError(sheet)
    a=c.abs();assert a.lo>0
    unit=c.conj()/a
    U=-2-F(2,3)*c*c
    V=F(8,3)*c-F(8,27)*c*c*c
    # Real-slice derivative data: E=unit e^(4i theta),
    # J=J0 e^(2i theta), C=C0 e^(i theta), Z=Z0 e^(2i theta).
    J=unit*U/2;B=unit*V/4;Z=U.conj()/2
    f=U.abs()/2;gg=V.abs()/4
    gamma=(I(3)/(8*a)).sqrt().sqrt()
    return (unit,J,B,Z),f,gg,gamma
