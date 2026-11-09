"""Degree-complete rational circle and signed exterior-coefficient checks for 9 through 14."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
import json
from coherent_projection_transfer import constants,rows
from audit_quintic_rational import I,PI
def run(n, nonlinear_bound):
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
            if n==9:extra=max(extra,(I0.l-1)*max(Q(0),b*c-nonlinear_bound*hi)**2)
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
        vals=run(n,Q(res['nonlinear_bound']['R9_upper']));res[str(n)]=vals;print(n,len(vals),'complete diameter intervals',flush=True)
    record = Path(__file__).with_name('joint_circle_intervals.json')
    record.write_text(json.dumps(res, separators=(',', ':')) + '\n',
                      encoding='utf-8', newline='\n')

