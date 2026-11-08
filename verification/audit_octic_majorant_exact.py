"""Elementary rational inequalities replacing the degree-eight unit-ball cover."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json

def verify():
    # After maximizing in v, the majorant is <= A(t)+B(t)w+k t^2 sqrt(w),
    # where t=sqrt(u), k=5sqrt(3)/16, and k^2=75/256.
    # A(t)<=a(t), B(t)<=1/40-5t/16.
    a=[F(5,14),F(5,16),F(5,896),-F(65,128),F(241,7168)]
    assert F(7,5)**2<2<F(10,7)**2 and 3<F(7,4)**2
    # 0<=t<=1/4: discard the negative cubic, use w<=1.
    small=sum(a[j]*F(1,4)**j for j in (0,1,2,4))+F(1,40)+F(35,64)*F(1,16)
    assert small<F(53,100)
    # 1/4<=t<=1/2: 25t-2>=17t, so Young's maximum costs <=375t^3/1088.
    cubic=F(65,128)-F(375,1088);assert cubic==F(355,2176)
    assert F(5,16)-3*cubic*F(1,4)>0
    middle=a[0]+a[1]/2+a[2]/4-cubic/8+a[4]/16
    assert middle<F(53,100)
    # 1/2<=t<=1: 25t-2>=21t. Bound t^2 and t^4 by t and t^3.
    linear=a[1]+a[2];negative=F(65,128)-F(375,1344)-a[4]
    assert linear==F(285,896) and negative==F(1399,7168)
    assert linear/(3*negative)==F(760,1399)<F(9,16)
    large=a[0]+F(2,3)*linear*F(3,4)
    assert large==F(925,1792)<F(53,100)
    assert max(small,middle)<large
    out={'status':'elementary rational bounds passed',
         'small_bound':str(small),'middle_bound':str(middle),
         'large_bound':str(large),'majorant_upper':'925/1792',
         'required_upper':'53/100','parameter_boxes':0}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Degree-eight majorant <=925/1792 <53/100; no adaptive boxes.')
    return out

if __name__=='__main__':assert __debug__;verify()
