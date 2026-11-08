"""Exact scalar algebra for the nonlinear ninth coefficient: |R9|<15/16.

The proof shares the quadratic budget with the cubic and quartic terms.
One univariate cubic remains; a quadratic lower bound proves its positivity.
No adaptive parameter cover or floating arithmetic is used.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json

def verify():
    # a=sqrt(12)/8 < 7/16; 5sqrt(2)/24 <5/16;
    # max (11sqrt(2)/32)(1-t)^(3/2)sqrt(t)=33sqrt(6)/512 <1/6.
    assert F(12,64)<F(7,16)**2
    assert 2<F(3,2)**2
    assert F(33,512)**2*6<F(1,6)**2
    # Completing the ellipse yields the bound
    # 7/16 +21s/160+s sqrt(623/800-45s/64).
    assert F(9,5)*(F(25,64)+F(121,576*5))==F(623,800)
    assert F(9,5)*F(25,64)==F(45,64)
    assert -F(7,32)+F(1,6)+F(2,5)*F(11,24)==F(21,160)
    p=[F(1,4),-F(21,160),-F(19495,25600),F(45,64)]
    def value(s):return sum(q*s**i for i,q in enumerate(p))
    def derivative(s):return sum(i*q*s**(i-1) for i,q in enumerate(p) if i)
    # p' < 0 for 0<=s<=1/2: bound its quadratic term by a linear term.
    assert 2*p[2]+3*p[3]/2==-F(1199,2560)
    assert value(F(1,2))==F(1677,20480)>0
    # On [1/2,1], p'' >=1501/2560. Use Taylor's theorem at 3/4.
    second_lower=2*p[2]+6*p[3]/2
    assert second_lower==F(1501,2560)>0
    assert value(F(3,4))==F(325,16384)
    assert derivative(F(3,4))==-F(891,10240)
    lower=value(F(3,4))-derivative(F(3,4))**2/(2*second_lower)
    assert lower==F(411311,30740480)>0
    assert F(1,2)-F(21,160)>0
    out={'status':'all elementary rational inequalities passed',
         'R9_upper':'15/16','positivity_lower':'411311/30740480',
         'parameter_subdivisions':0}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Nonic nonlinear coefficient: |R9|<15/16 by a quadratic lower bound.')
    return out

if __name__=='__main__':assert __debug__;verify()
