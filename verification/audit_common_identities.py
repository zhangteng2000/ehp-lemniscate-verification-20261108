"""Exact checks of algebraic identities to be used in the proposed refactoring.

These checks audit identities, not the EHP comparison or a parameter cover.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction
import json
import sympy as sp

X,Y,U,V,c,h=sp.symbols('X Y U V c h', real=True)
A,B,C,D=Y,3*h*X,V-3*h*h*Y,c+h*U-h**3*X
discriminant=18*A*B*C*D-4*B**3*D+B**2*C**2-4*A*C**3-27*A**2*D**2
zr,zi=X*U+Y*V,Y*U-X*V
proposed=(h*zi+Y*c)**2+sp.Rational(4,27)*Y*V**3+h*h*(4*h*X*c-sp.Rational(4,3)*V*V+4*h*h*(zr-h*h))
rem=sp.rem(sp.Poly(sp.expand(discriminant+27*proposed),X),sp.Poly(X*X+Y*Y-1,X))
assert rem.is_zero

# Check the coefficient formula for every pair J,l and m<=8.
# The shape of the formula is independent of n; x enters only through b.
t,u=sp.symbols('t u')
def coefficients(f,N):
    return sp.Poly(sp.series(f,t,0,N+1).removeO(),t)

checks=0
for b in (sp.Rational(2,5),sp.Rational(5,12),sp.Rational(3,7)):
    for J in range(1,5):
        for l in range(0,4):
            N=8-J-l
            if N<0:continue
            # The part left after t^(J+l)(1-u)^J u^l is removed.
            first=coefficients((1-u*t)**(b-J-l),N)
            second=coefficients((1-t)**(-b-J),N)
            for r in range(N+1):
                actual=sp.expand(sum(first.nth(k)*second.nth(r-k) for k in range(r+1)))
                expected=sum(sp.rf(J+l-b,k)/sp.factorial(k)*sp.rf(b+J,r-k)/sp.factorial(r-k)*u**k for k in range(r+1))
                assert sp.expand(actual-expected)==0
                checks+=1
    N=8
    first=coefficients((1-u*t)**b,N)
    second=coefficients((1-t)**(-b),N)
    for m in range(1,N+1):
        numer=sp.expand(sum(first.nth(k)*second.nth(m-k) for k in range(m+1)))
        assert sp.rem(sp.Poly(numer,u),sp.Poly(1-u,u)).is_zero

# Half-integer beta coefficients in the two-weight energy proof.
for j in range(9):
    for l in range(9-j):
        p=j+l
        def beta(a,b):return sp.gamma(a)*sp.gamma(b)/sp.gamma(a+b)
        def cc(j,l):
            p=j+l
            return sp.binomial(2*j,j)/(j+1)/(2*sp.pi*4**j)*sp.binomial(2*p,2*l)*beta(l+sp.Rational(1,2),j+sp.Rational(3,2))*beta(sp.Rational(1,2),p+sp.Rational(1,2))
        def bb(k):return sp.binomial(2*k,k)/4**k
        assert sp.simplify(cc(j,l)-sp.binomial(p,j)*cc(p,0)*bb(j+1)/bb(p+1))==0

out={'cubic_discriminant':'exact polynomial identity modulo X^2+Y^2=1',
     'exterior_coefficient_checks':checks,
     'model_numerator_divisibility':'all m=1,...,8, three rational exponents',
     'pair_beta_coefficient_checks':45,
     'scope':'algebraic checks only; does not certify any parameter coverage'}
Path(__file__).with_name('common_identities_audit.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
