"""Exact symbolic extraction of the two Laurent coefficients used in the paper.

This checks the finite affine-flow algebra. The analytic stationarity,
boundary Parseval, and orthogonality arguments are proved in the manuscript.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

import sympy as S

def verify():
    b,u,c,d,h,eta,h2=S.symbols('b u c d h eta h2')
    sb=b*(b+1)/2;tb=b*b*(b-1)*(b+3)/4
    p=c*c*(sb+b*(1-b)*u/2)+d*c*h*(1+b+(1-b)*u)+d*d*eta*u+d*d*h2*(1-u)
    series=1-b*u+b*(b-1)*u*u/2
    got=S.expand(series*p).coeff(u,2)
    expected=tb*(c*c+2*d*c*h/b)-b*d*d*eta+sb*d*d*h2
    assert S.factor(got-expected)==0
    # The first Laurent row comes from (1-u)^b h_n(p).
    assert S.diff((1-u)**b,u).subs(u,0)==-b
    print('Exact affine-flow second coefficient and first-row constant: passed.')

if __name__=='__main__':verify()
