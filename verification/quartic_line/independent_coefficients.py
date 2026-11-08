"""Independent rational coefficient enclosures for the four quartic charts.

The first generator rationalizes the imaginary radical and rounds each
operation to a dyadic interval.  This generator keeps all rational operations
exact, takes square roots on a 100-bit grid, and uses the monotonicity of
u(s)^2=2 sqrt(1+s+s^2)+2+s and v(s)^2=2 sqrt(1+s+s^2)-2-s.
No floating-point operation contributes to an acceptance test.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json

S = 1 << 100

class R:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        assert self.lo <= self.hi
    def __add__(self, b):
        b = interval(b); return R(self.lo+b.lo, self.hi+b.hi)
    __radd__ = __add__
    def __neg__(self): return R(-self.hi, -self.lo)
    def __sub__(self, b): return self + -interval(b)
    def __rsub__(self, b): return interval(b) + -self
    def __mul__(self, b):
        b = interval(b)
        v = [x*y for x in (self.lo,self.hi) for y in (b.lo,b.hi)]
        return R(min(v),max(v))
    __rmul__ = __mul__
    def __truediv__(self, b):
        b=interval(b); assert b.lo*b.hi > 0
        return self*R(1/b.hi,1/b.lo)
    def __rtruediv__(self, b): return interval(b)/self
    def square(self):
        v=(self.lo**2,self.hi**2)
        return R(0 if self.lo <= 0 <= self.hi else min(v),max(v))
    def sqrt(self):
        assert self.lo >= 0
        a=isqrt(self.lo.numerator*S*S//self.lo.denominator)
        b=isqrt(self.hi.numerator*S*S//self.hi.denominator)
        assert F(a*a,S*S) <= self.lo < F((a+1)**2,S*S)
        assert F(b*b,S*S) <= self.hi < F((b+1)**2,S*S)
        return R(F(a,S),F(b if F(b*b,S*S)==self.hi else b+1,S))
    def dyadic(self,bits):
        return (self.lo*(1<<bits)).__floor__(),(self.hi*(1<<bits)).__ceil__()

def interval(x): return x if isinstance(x,R) else R(x)

class Z:
    __slots__=('x','y')
    def __init__(self,x=0,y=0):self.x,self.y=interval(x),interval(y)
    def __add__(self,b):
        b=complex_interval(b);return Z(self.x+b.x,self.y+b.y)
    __radd__=__add__
    def __neg__(self):return Z(-self.x,-self.y)
    def __sub__(self,b):return self+-complex_interval(b)
    def __rsub__(self,b):return complex_interval(b)+-self
    def __mul__(self,b):
        b=complex_interval(b)
        return Z(self.x*b.x-self.y*b.y,self.x*b.y+self.y*b.x)
    __rmul__=__mul__
    def __truediv__(self,b):return Z(self.x/b,self.y/b)
    def bar(self):return Z(self.x,-self.y)
    def norm(self):return (self.x.square()+self.y.square()).sqrt()

def complex_interval(x):return x if isinstance(x,Z) else Z(x)

def radicals(s):
    r=(R(1)+s+s*s).sqrt()
    u=(2*r+2+s).sqrt()
    w=2*r-2-s
    # v(0)=0 exactly; outward subtraction may give a negative lower endpoint.
    v=R(max(F(0),w.lo),max(F(0),w.hi)).sqrt()
    return u,v

def build(parts=128,bits=24,output_dir=None):
    target=Path(__file__).resolve().parent if output_dir is None else Path(output_dir)
    rows=[];meta=[]
    for sheet in range(3):
        for j in range(parts):
            l,h=F(j,parts),F(j+1,parts)
            s0,s1=1-h*h,1-l*l
            u0,v0=radicals(s0);u1,v1=radicals(s1)
            u=R(u0.lo,u1.hi);v=R(v0.lo,v1.hi);t=R(l,h)
            if sheet==0:c=Z(t+u)
            elif sheet==1:c=Z(t-u)
            elif sheet==2:c=Z(-t,v)
            else:c=Z(-t,-v)
            a=c.norm();assert a.lo>0
            e=c.bar()/a
            U=-2-F(2,3)*c*c
            V=F(8,3)*c-F(8,27)*c*c*c
            data=(e,e*U/2,e*V/4,U.bar()/2)
            nums=[]
            for z in data:
                for q in (z.x,z.y):nums.extend(q.dyadic(bits))
            nums.extend(((U.norm()/2).dyadic(bits)[1],(V.norm()/4).dyadic(bits)[1]))
            gamma=(R(3)/(8*a)).sqrt().sqrt()
            ident=sheet*parts+j
            rows.append(' '.join(map(str,[ident,*nums])))
            meta.append({'id':ident,'sheet':sheet,'left':str(l),'right':str(h),
                         'gamma_lower':str(gamma.lo)})
    path=target/f'independent_coefficients_{parts}.txt'
    path.write_text('\n'.join(rows)+'\n',encoding='ascii')
    (target/f'independent_parameters_{parts}.json').write_text(json.dumps(
        {'intervals_per_sheet':parts,'coefficient_bits':bits,'rows':meta},indent=1)+'\n',encoding='utf-8')
    print('Independently enclosed',len(rows),'closed parameter intervals.')
    return path

if __name__=='__main__':
    assert __debug__;build()
