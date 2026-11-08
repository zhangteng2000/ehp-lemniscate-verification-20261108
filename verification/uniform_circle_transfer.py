"""Exact endpoint arithmetic for the analytic circle implication in every n>=28.

There is no degree loop and no parameter cover. The monotonicity argument
is analytic; the inequalities checked below are the rational endpoint
inequalities explicitly displayed in the manuscript.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json

def atan_bounds(q,N):
    assert N%2==0
    lo=sum(F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N))
    return lo,lo+F(1,(2*N+1)*q**(2*N+1))

def verify():
    l,u=atan_bounds(5,60);a,b=atan_bounds(239,16)
    plo,phi=F(333,106),F(355,113)
    assert plo<16*l-4*b<16*u-4*a<phi
    ln2lo=2*sum(F(1,(2*j+1)*3**(2*j+1)) for j in range(12))
    assert 2*phi-8*ln2lo<F(739,1000)
    assert 7*ln2lo>2  # monotonicity of (128 n^2)^(-1/n).
    q=F(739,800);bs=F(27,56);a=F(543,5000);z=F(299,100)
    assert q**28<a and F(28)/plo<z*z
    A=2*plo*(bs*bs+(bs*(1-bs)/2)**2)
    assert A>F(31,20)
    assert 1-a*a/16>F(999,1000)**2
    assert 1-a*a>F(497,500)**2
    assert F(33,50)**28<F(1,128*28*28)
    assert F(19,20)**28<F(1,4)
    assert bs*F(999,1000)-z*a/4>F(2,5)
    assert bs*F(497,500)-z*a>F(309,2000)
    assert (1-q)*56>1  # monotonicity of sqrt(n) q^n.
    first=F(31,20)*F(2,5)**2+F(4,5)*F(33,50)
    second=F(31,20)*F(309,2000)**2+F(739,1000)*F(19,20)
    assert first==F(97,125)>F(739,1000)
    assert second-F(739,1000)==F(3911,80000000)>0
    assert F(80,3)**6<128**5  # 1/2+8/128^(5/6)<4/5.
    out={'status':'all uniform analytic endpoint inequalities passed',
         'minimum_degree':28,'degree_subdivisions':0,'diameter_partitions':0,
         'first_case_lower':'97/125','second_case_margin':'3911/80000000'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Uniform circle implication for all n>=28: endpoint arithmetic passed; no degree or diameter cover.')
    return out

if __name__=='__main__':verify()
