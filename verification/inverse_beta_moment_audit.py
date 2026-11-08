"""Exact arithmetic in the centered inverse-beta lemma.

The universal recurrence and its integration-by-parts proof are in the
manuscript. The raw-moment comparisons below are extra identity checks;
the acceptance bounds for the infinite degree range are the rational
endpoint and monotonicity inequalities, not samples of degrees.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

def verify():
    u=F(33,256);W=F(46,5);m0=223
    A=2*(1+u)*W/m0;B=2*u*W*W/(3*m0)
    assert A/4+B<F(1,16)
    low=2*u*W*W/(3*222)
    assert low<F(1,30)
    assert W/((1-u)*(F(m0,4)-F(1,2)))<1
    assert 12*W/((1-u)**2*256)<F(3,5)
    assert 9*W/m0<F(2,5)
    assert 1-u>F(3,4) and 2**4<3**3
    # Differentiating n r^(c n): log r <= -(1-r) makes both
    # logarithmic derivatives negative for all n>=256.
    assert F(1,256)<(1-u)/5
    assert F(1,256)<3*(1-u)/5
    tail=256*(F(3,5)**111+F(2,5)**223)
    assert tail<F(1,1000)
    assert low+tail<F(1,25)
    checks=0
    for n,k in ((256,2),(256,16),(256,32),(256,33),(257,33),(512,46),(1024,64)):
        m=n-k;N=n-1;a=F(N,m);p=F(k-1,N)
        assert 2<=k<=u*n and m>=m0
        moments=[F(1),F(0)]
        for j in range(2,25):
            moments.append(F(j-1,m-j+1)*(p*moments[-2]+(1+p)*moments[-1]))
            raw=sum(F(comb(j,r))*(-1)**(j-r)*a**(-r)*F(comb(N,r),comb(m,r)) for r in range(j+1))
            assert raw==moments[-1]>=0
            checks+=1
    out={'status':'all displayed endpoint inequalities and supplemental raw-moment identities passed',
         'low_order_sum_upper':str(low),'tail_sum_lt':'1/1000',
         'total_sum_lt':'1/25','minimum_degree':256,
         'raw_moment_identity_checks':checks,'degree_sampling_used_for_uniform_bound':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Centered inverse-beta lemma: all rational endpoint bounds passed; recurrence identities',checks)
    return out

if __name__=='__main__':verify()
