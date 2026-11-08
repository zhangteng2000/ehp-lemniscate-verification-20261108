"""Independent rational validation of all radius witnesses in the seeded diameter-two reduction."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import json
from audit_conditional_line import band,cond
from seeded_coefficient_descent import descend
G=Path(__file__).resolve().parent
@lru_cache(None)
def setup(n,k):
    m=n-k;pars=json.loads((G/'conditional_contours.json').read_text());out=[]
    for (l,h,v),tc in zip(band(m),pars[str(n)][str(k)],strict=True):
        K,E=cond(n,k,l,v,Q(tc,100));out.append((l,h,v,K,E))
    return out
def validate(n,L,start,lo,record):
    c=start.copy();assert len(record)==n-3
    for k,data in zip(range(n-2,1,-1),record,strict=True):
        m=n-k;ratio=[Q(comb(n,j),comb(m,j)) for j in range(m+1)]
        trial=[];witness=list(map(Q,data));arr=setup(n,k);assert len(witness)==len(arr)
        for A,(a,b,e,K,E) in zip(witness,arr):
            assert A>0
            if (K-b*E)*A<L:
                assert a>0
                r=a*A;sum_co=Q(0)
                for j in range(2,m+1):
                    bound=e[j]*A**j
                    if c[n-j] is not None:bound=min(bound,c[n-j]/ratio[j])
                    sum_co+=bound*r**(m-j)
                assert sum_co<r**m,(n,k,a,A)
            elt=[ratio[j]*e[j]*A**j for j in range(2,m+1)]
            if lo is None or all(lo[n-j]<=t for j,t in zip(range(2,m+1),elt)):
                trial.append(elt)
        if not trial:
            assert lo is not None and any(lo[n-j]>0 for j in range(2,m+1))
            return None
        for j,part in enumerate(zip(*trial),2):
            t=max(part);c[n-j]=t if c[n-j] is None else min(t,c[n-j])
    return c
def passes(n,L,start,lower,rounds=2):
    # The first program only suggests radii. The accepting checks use separate
    # envelopes, trig estimates, series products, contour constants and caps.
    pars=json.loads((G/'conditional_contours.json').read_text());c=start.copy()
    for _ in range(rounds):
        before=c.copy()
        witness=descend(n,L,pars,c,lower=lower)
        if witness is None:
            # Recover all bins including an impossible step with an explicit independent check.
            # The raw generator stopped partway; rerun its arithmetic while collecting the prefix.
            w,ended=prefix(n,L,pars,before,lower)
            assert ended
            checked=validate_partial(n,L,before,lower,w);assert checked is None
            return None
        independently=validate(n,L,before,lower,witness)
        if independently is None:return None
        assert all(independently[j]<=c[j] for j in range(2,n+1))
    return c
def prefix(n,L,pars,c,lower):
    # Collect proposed witnesses without making any acceptance decision here.
    from verify_conditional_line import make_envelope,errors,up
    records=[];c=c.copy()
    for k in range(n-2,1,-1):
        m=n-k;row=[];new=[]
        for (a,b,e),arg in zip(make_envelope(m,128),pars[str(n)][str(k)],strict=True):
            K,E=errors(n,k,(a,b,e),arg);g=up(L/(K-E*b)) if K>E*b else None
            full=[Q(comb(n,j),comb(m,j))*e[j] for j in range(2,m+1)]
            def fit(A):
                return a!=0 and sum((w if c[n-j] is None else min(w,c[n-j]/A**j))*Q(comb(m,j),comb(n,j))/a**j
                                   for j,w in enumerate(full,2))<1
            high=Q(32) if g is None else min(g,Q(32))
            if fit(high):
                low=0
                for z in range(18):
                    mid=(low+high)/2
                    if fit(mid):high=mid
                    else:low=mid
                g=up(high)
            assert g is not None
            row.append(g);vec=[up(w*g**j) for j,w in enumerate(full,2)]
            if all(lower[n-j]<=v for j,v in enumerate(vec,2)):new.append(vec)
        records.append(row)
        if not new:return records,True
        for j,v in enumerate(zip(*new),2):
            b=max(v);c[n-j]=b if c[n-j] is None else min(b,c[n-j])
    return records,False
def validate_partial(n,L,c,lower,ws):
    # The unknown slots after the impossible step are never reached.
    padded=ws+[[] for _ in range(n-3-len(ws))]
    return validate(n,L,c,lower,padded)

