"""A rational absolute-error enclosure of both binary64 circuits in QuinticPhaseAudit.cs."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json

U=F(1,2**52)
SMALL=F(1,2**300)

class V:
    # Absolute bound on the exact real number and on the computed error.
    def __init__(self,a,e=0):self.a,self.e=abs(F(a)),F(e)
    def __add__(self,o):
        o=cv(o);m=self.a+o.a;v=self.e+o.e
        return V(m,v+U*(m+v)+SMALL)
    __radd__=__add__
    __sub__=__add__
    __rsub__=__add__
    def __mul__(self,o):
        o=cv(o);m=self.a*o.a;v=self.a*o.e+o.a*self.e+self.e*o.e
        return V(m,v+U*(m+v)+SMALL)
    __rmul__=__mul__

def cv(a):return a if isinstance(a,V) else V(a)
def rnd(a):a=F(a);return V(a,a*U+SMALL)

f=V(F(9,4));g=V(1);h=V(F(7,4));h2=h*h
tr=V(1,F(1,2**46))
y,x=tr,tr;bi=f*tr;c=g*tr;zi=f*tr;zr=f*tr
yc=y*c;free=rnd(F(4,27))*y*bi*bi*bi
odd=4*x*c;even=rnd(F(4,3))*bi*bi
s=h*zi+yc
D=s*s+free+h2*(h*odd-even+4*h2*(zr-h2))

# Shape parameters, steps, |h| and both extrema of |h| are all exactly
# represented dyadic values on every grid used here. max(0,H) is Lipschitz.
hi=lo=h;r2=hi*hi;r3=r2*hi;r4=r2*r2;f2=f*f
w=g*g+rnd(F(4,27))*f2*f
ang=324*w+242*f*g*hi+64*f2*r2+324*g*r3+16*f*r4
alpha=16*w+8*f*g*hi+rnd(F(16,3))*f2*r2+16*g*r3
beta=4*g*g+2*f*g*hi+4*g*r3
shape=2*r2+rnd(F(8,9))*f
height=48*f*r2+24*g*hi+2*f2-120*(lo*lo)*(lo*lo)
dt,da,db,ds,dh=map(V,(F(1,256),F(1,32),F(1,32),F(9,32),F(7,256)))
er=(F(79,8)*(dt*dt*ang+da*da*alpha+db*db*beta)+shape*ds*ds+height*dh*dh)*F(1,8)
withpad=er+F(1,2**25)
total=D.e+withpad.e
assert total<F(1,2**28)
print("Absolute combined roundoff and trig error <",float(total*F(1000001,1000000)),"< 2^-28 < guard 2^-25")
root=Path(__file__).resolve().parent
(root/"sextic_roundoff_audit.json").write_text(json.dumps({"combined_error_lt":"2^-28","guard":"2^-25"},indent=2)+"\n")

