"""Rational upper bound retaining the favorable low exterior diagonals."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from audit_quintic_rational import I,PI,J

def lower_dyadic(z,bits=28):return F((z*(1<<bits)).__floor__(),1<<bits)
def upper_dyadic(z,bits=28):return F((z*(1<<bits)).__ceil__(),1<<bits)

def certified(n,R,t,b,constants,hn_bound=None,components=False):
    x,_,C,_,U,_=constants
    H=(U/(2*PI)-1)**F(1,2)
    R=F(R);t=F(t);b=F(b);rc=1/(1-b);u=t*R
    assert rc<R and t<b and u<=F(1,2)
    S=I(R)/(1+I(u*u))**F(1,2)
    s=lower_dyadic(S.l); y=upper_dyadic((S**(-x)).h)
    assert s>1 and y<1
    eps=u/(1-I(u*u))**F(1,2)
    core=U*(I(R)/2)**x*2*J(1-x)/PI*(eps+eps**2)+(1+x)/(2+x)*t*t*I(R)**(2+x)
    r=(rc/R)**4
    lin=H*I(R)**x*(rc/R)**2/((2-x)*(1-I(r))**F(1,2))
    quad=2*H**2*I(rc)**(4*x)*I(R)**(-3*x)/(3*x*(1-r))
    move=(4*lin+2*quad)*t/(b-t)
    hn=H.h
    if n==5:hn=min(hn,(F(1,5)+1/(4*I(2)**F(1,2))).h+t/5)
    if n==6:hn=min(hn,F(21,40)+t/6)
    if hn_bound is not None:hn=min(hn,F(hn_bound))
    first=4*hn*J(x)*I(t)**(1-x)
    top=16*n+64
    def K(p,l):
        return 1+1/p+(1+F(11,7)*l)/min(p-1,1)
    powers=[y**i for i in range(2*top+1)]
    sums=[]
    for i in range(2,top):
        sm=F(0)
        for j in range(i%n,top,n):
            if j<2 or i+j<=n+1:continue
            p=F(i+j-1,n);l=abs(i-j)//n
            sm+=K(p,l)*powers[i+j]
        sums.append(sm)
    # Absolute sum of EVERY discarded matrix coefficient bounds each discarded row.
    tail=2*y**(top+2)*((n+6+2*top)/(1-y)**2+4*y/(1-y)**3)
    row=max(sums)+tail
    V=F(1,s);Q=34
    linear=sum((K(F(q)-x,q)*V**q)**2 for q in range(2,Q))
    extra=V**Q*((3+2*Q)/(1-V)+2*V/(1-V)**2)
    lin2=I(linear)**F(1,2)+extra
    high=2*t*I(R)**(1+x)*(H**2*row+2*H*lin2)
    if components:
        return {'core':core.h,'motion':move.h,'first':first.h,'higher':high.h,
                'total':(core+move+first+high).h,'radius_lower':s,'root_bound':hn,
                'matrix_row_upper':row,'linear_tail_upper':lin2.h}
    return (core+move+first+high).h

