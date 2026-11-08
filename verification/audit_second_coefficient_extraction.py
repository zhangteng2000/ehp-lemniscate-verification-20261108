"""Supplemental exact checks for the shorter second-coefficient derivation.

The universal coefficient extraction is proved analytically in the paper.
These finite direct polynomial-flow checks do not replace that proof and
are not inputs to any numerical acceptance test.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
from audit_stationary_direct import (Dual, add, mul, derivative, divide,
                                     inverse, constant_after_inverse)

def verify():
    count=0
    for n in range(4,10):
        b=Q(n-1,2*n);sb=b*(b+1)/2
        for seed in (1,2):
            q=[Q((-1)**(j+seed)*(j+seed),n+2) for j in range(n-1)]+[Q(0),Q(1)]
            _,r=divide(mul(q,q),derivative(q))
            _,hq=inverse([Dual(v,r[j] if j<len(r) else 0) for j,v in enumerate(q)],2*n)
            c=Q(seed,3);d=Q(1,2);delta=d**n
            p=[delta*q[j]/d**j for j in range(n+1)];p[0]+=c
            z,w=divide(mul(p,p),derivative(p))
            zdot,_=divide(add([2*v for v in mul(p,w)],
                              [-v for v in mul(derivative(w),z)]),derivative(p))
            a2=add(add([v/4 for v in derivative(zdot)],
                       [v/4 for v in mul(derivative(derivative(z)),z)]),
                   [v/8 for v in mul(derivative(z),derivative(z))])
            ci,hp=inverse(p,2*n)
            _,variation=inverse([Dual(v,w[j] if j<len(w) else 0) for j,v in enumerate(p)],n)
            eta_p=variation[n].b
            gamma=constant_after_inverse(a2,ci,hp)
            assert gamma==sb*hp[2*n]-b*eta_p
            assert hp[2*n]==sb*c*c+(b+1)*delta*c*hq[n].a+delta*delta*hq[2*n].a
            assert eta_p==b*c*c+2*delta*c*hq[n].a+delta*delta*hq[n].b
            count+=1
    here=Path(__file__).resolve()
    receipt={'checked_utc':datetime.now(timezone.utc).isoformat(),
             'cases':count,'exact_equalities':3*count,'exit_code':0,
             'source_sha256':sha256(here.read_bytes()).hexdigest(),
             'primitive_sha256':sha256(here.with_name('audit_stationary_direct.py').read_bytes()).hexdigest(),
             'scope':'Supplemental exact rational examples for the analytic universal and affine coefficient identities; not a proof for arbitrary degree.'}
    here.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('Universal and affine second-coefficient extraction:',count,'exact cases,',3*count,'equalities passed.')

if __name__=='__main__':verify()
