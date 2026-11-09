"""Degree-complete rational circle and signed exterior-coefficient checks for 9 through 14."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import json
from coherent_projection_transfer import constants,rows
from audit_quintic_rational import I,PI
def nonlinear():
    # 9/2-99(t+t^2)/70 is positive. Its square minus 16(1+4t-5t^2)/3.
    S=[Q(9,2),-Q(99,70),-Q(99,70)];p=[Q(0)]*5
    for i,a in enumerate(S):
        for j,b in enumerate(S):p[i+j]+=a*b
    for j,v in enumerate((1,4,-5)):p[j]-=Q(16,3)*v
    o=[]
    for i in range(4):
        z=Q(i,4);w=Q(1,4)
        co=[sum(p[j]*comb(j,l)*z**(j-l)*w**l for j in range(l,5)) for l in range(5)]
        be=[sum(co[j]*Q(comb(l,j),comb(4,j)) for j in range(l+1)) for l in range(5)]
        assert min(be)>0;o.append(min(be))
    assert Q(99,70)**2>2
    # quadratic: sqrt(12)/8, cubic: (5/48)(9/2), quartic: 33sqrt(6)/512.
    assert ((I(12)**Q(1,2))/8+Q(5,48)*Q(9,2)+33*I(6)**Q(1,2)/512).h<Q(17,16)
    return list(map(str,o))
def run(n):
    prop=constants(n);x,b,C,I0,U,_=prop
    L=Q(41,50) if n==9 else Q(7,8)
    data=rows(n,9 if n==9 else 4,prop)
    lo0=Q(1,128*n*n);stack=[(lo0,Q(1,2)),(Q(1,2),Q(1))];ans=[]
    while stack:
        lo,hi=stack.pop();extra=Q(0)
        if hi<=Q(1,2):
            c=(I(1-hi**2)**Q(1,2)).l;cb=(I(1-lo**2)**Q(1,2)).h
            for j,(N,band) in enumerate(data,1):
                sub=sum(s*hi**i*cb**(j-i) for i,s in band.items())
                extra+=max(Q(0),c**j*N.l-sub)**2/N.h
            if n==9:extra=max(extra,(I0.l-1)*max(Q(0),b*c-Q(17,16)*hi)**2)
        target=((U-C).h-2*PI.l*extra)/L
        if target<=0 or target**n<lo:
            ans.append((lo,hi))
        else:
            assert hi-lo>Q(1,2**28),(n,float(lo),float(hi))
            mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
    ans.sort()
    assert ans[0][0]==lo0 and ans[-1][1]==1 and all(a[1]==b[0] for a,b in zip(ans,ans[1:]))
    return [[str(a),str(b)] for a,b in ans]
if __name__=='__main__':
    from audit_nonic_majorant_exact import verify as verify_nonlinear
    assert __debug__;res={'nonlinear_bound':verify_nonlinear()}
    for n in range(9,15):
        vals=run(n);res[str(n)]=vals;print(n,len(vals),'complete diameter intervals',flush=True)
    record = Path(__file__).with_name('joint_circle_intervals.json')
    record.write_text(json.dumps(res, separators=(',', ':')) + '\n',
                      encoding='utf-8', newline='\n')

