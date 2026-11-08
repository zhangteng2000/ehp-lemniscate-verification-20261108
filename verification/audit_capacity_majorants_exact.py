"""Elementary quadratic maxima replace two adaptive box covers.

The analytic reduction to these polynomials is in the manuscript. This file
checks only the displayed finite algebra; it is not a substitute for that
reduction. No root approximation, external package, or floating arithmetic is
used for acceptance.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json

def verify():
    # Degree six: the two ranges u<=4/7 and u>=4/7.
    assert F(4,7)<F(4,5)**2
    assert F(1,4)+F(4,7)/48<F(17,64)<F(27,100)
    # u^(3/2)<=(5u-1)/4 on u>=4/7 follows after t=sqrt(u):
    # (5t^2-1)-4t^3 = (1-t)(4t^2-t-1), whose last factor is positive.
    assert 4*F(4,7)-1-1>0  # t<=1, so 4t^2-t-1 >= 16/7-2>0.
    assert 7537<87**2 and F(15+87,384)==F(17,64)
    # Degree seven: the quartic is 9/25 minus an explicit upper bound.
    assert F(7,5)**2<2<F(10,7)**2
    power=[F(11,100),-F(25,112),-F(57,256),F(35,64),-F(25,512)]
    def val(t):return sum(q*t**i for i,q in enumerate(power))
    def dval(t):return sum(i*q*t**(i-1) for i,q in enumerate(power) if i)
    a=F(1,4)
    coarse=F(11,100)-F(25,112)*a-F(57,256)*a*a-F(25,512)*a**4
    assert coarse>F(1,25)
    second=2*power[2]+6*power[3]*a+12*power[4]*a*a
    assert second==F(693,2048)>0
    assert 2*power[2]+6*power[3]+12*power[4]>second
    # q'' is concave, so these endpoint bounds apply on [1/4,1].
    t=F(9,16);lower=val(t)-dval(t)**2/(2*second)
    assert lower==F(1428481870703,227881569484800)>F(3,500)
    out={'status':'all displayed exact scalar inequalities passed',
         'degree_six_bound':'17/64 < 27/100',
         'degree_seven_bound':'9/25','quartic_lower':str(lower),
         'adaptive_boxes':0}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Capacity majorants: degree six <=17/64; degree seven <9/25, exact algebra passed.')
    return out

if __name__=='__main__':assert __debug__;verify()
