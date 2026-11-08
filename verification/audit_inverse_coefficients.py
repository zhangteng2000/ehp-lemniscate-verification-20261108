"""Exact algebra for the shared logarithmic identity and its small-degree terms.

These symbolic checks corroborate, but do not replace, the analytic proof of
the inverse germ, its area bound, and the target-curve comparison.
"""
from pathlib import Path
import json
import sympy as s

if not __debug__:
    raise RuntimeError('Proof verification requires Python assertions; do not use -O.')

def verify():
    t=s.symbols('t');a=s.symbols('a1:9')
    Z=sum(j*a[j-1]*t**(j+1) for j in range(1,9))
    R=s.series(s.sqrt(1-Z)-1-s.log(1-Z)/2,t,0,10).removeO().expand()
    expected={
        3:0,
        4:a[0]**2/8,
        5:a[0]*a[1]/2,
        6:s.Rational(3,4)*a[0]*a[2]+a[1]**2/2+s.Rational(5,48)*a[0]**3,
        7:a[0]*a[3]+s.Rational(3,2)*a[1]*a[2]+s.Rational(5,8)*a[0]**2*a[1],
        8:s.Rational(5,4)*a[0]*a[4]+2*a[1]*a[3]+s.Rational(9,8)*a[2]**2
          +s.Rational(15,16)*a[0]**2*a[2]+s.Rational(5,4)*a[0]*a[1]**2+s.Rational(11,128)*a[0]**4,
        9:s.Rational(3,2)*a[0]*a[5]+s.Rational(5,2)*a[1]*a[4]+3*a[2]*a[3]
          +s.Rational(5,4)*a[0]**2*a[3]+s.Rational(15,4)*a[0]*a[1]*a[2]
          +s.Rational(5,6)*a[1]**3+s.Rational(11,16)*a[0]**3*a[1]
    }
    for n,want in expected.items():
        assert s.expand(R.coeff(t,n)-want)==0
        for powers,co in s.Poly(R.coeff(t,n),*a).terms():
            if co==0:continue
            assert co>=0 and sum(powers)>=2
            assert sum((j+2)*v for j,v in enumerate(powers))==n
    r,tau=s.symbols('r tau',real=True,nonzero=True)
    M=lambda z:z/(1-s.I*tau*z)
    assert s.simplify(1/M(r)+1/M(-r)+2*s.I*tau)==0
    # Check the resultant sign in both parities without a square-root convention.
    z,w=s.symbols('z w')
    for n in range(3,10):
        P=z**n+2*z**2+z-3
        lhs=s.resultant(P-w,s.diff(P,z),z)
        # Swapping degrees n and n-1 has even sign; the explicit (-1)^(n-1)
        # here instead distinguishes products (v_j-w) from (w-v_j).
        swapped=s.resultant(s.diff(P,z),P-w,z)
        assert s.expand(lhs-swapped)==0
        vals=s.Poly(swapped,w)
        assert vals.LC()==(-1)**(n-1)*n**n
    out={'status':'all exact inverse-coefficient identities passed',
         'degrees':list(expected),'signed_tail_identity':'1/M(r)+1/M(-r)=-2i tau',
         'resultant_parities_checked':list(range(3,10)),
         'scope':'symbolic algebra; analytic domains and coverage are checked separately'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Exact R3,...,R9, homogeneity, resultant signs, and signed tail identity passed.')
    return out

if __name__=='__main__':verify()
