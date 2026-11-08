"""Forward the cubic strip's source bounds into EVERY lower derivative's root-radius test."""
from verify_conditional_line import make_envelope,errors,choose,up,Q
from math import comb
from pathlib import Path
import json,argparse
G=Path(__file__).resolve().parent
def descend(n,L,params,co,parts=128,lower=None):
    rec=[]
    for k in range(n-2,1,-1):
        m=n-k;en=make_envelope(m,parts);fac=[Q(comb(n,j),comb(m,j)) for j in range(m+1)]
        cur=[None if co[n-j] is None else co[n-j]/fac[j] for j in range(m+1)]
        new=[Q(0)]*(m+1);ws=[];viable=0
        for (a,b,v),T in zip(en,params[str(n)][str(k)],strict=True):
            K,E=errors(n,k,(a,b,v),T);ub=up(L/(K-b*E)) if K>b*E else None
            def good(B):
                return a>0 and sum(min(v[j],cur[j]/B**j) / a**j if cur[j] is not None else v[j]/a**j
                                   for j in range(2,m+1))<1
            hi=Q(32) if ub is None else min(ub,Q(32))
            if good(hi):
                lo=Q(0)
                for _ in range(18):
                    mid=(lo+hi)/2
                    if good(mid):hi=mid
                    else:lo=mid
                ub=up(hi)
            assert ub and (ub*(K-b*E)>=L or good(ub))
            ws.append(str(ub))
            candidate=[up(v[j]*ub**j*fac[j]) for j in range(2,m+1)]
            if lower is not None and any(lower[n-j]>bound for j,bound in zip(range(2,m+1),candidate)):continue
            viable+=1
            for j,bound in zip(range(2,m+1),candidate):new[j]=max(new[j],bound)
        rec.append(ws)
        if lower is not None and not viable:
            assert any(lower[n-j]>0 for j in range(2,m+1));return None
        for j in range(2,m+1):co[n-j]=new[j] if co[n-j] is None else min(co[n-j],new[j])
    return rec
def run(n,L,f,g,nit=1):
    co=[None]*(n-1)+[Q(0),Q(1)]
    co[n-2]=Q(comb(n,2),3)*(L/f)**2;co[n-3]=comb(n,3)*(L/g)**3
    pars=json.loads((G/'conditional_contours.json').read_text())
    for _ in range(nit):ws=descend(n,L,pars,co)
    return {'coefficients':[None if c is None else str(c) for c in co],'energy':ws}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--L',type=Q,required=True)
    p.add_argument('--f',type=Q,required=True);p.add_argument('--g',type=Q,required=True);p.add_argument('--rounds',type=int,default=2)
    p.add_argument('--out',required=True);a=p.parse_args();r=run(a.n,a.L,a.f,a.g,a.rounds)
    print(a.n,a.f,a.g,[None if x is None else round(float(Q(x)),4) for x in r['coefficients']])
    (G/a.out).write_text(json.dumps({str(a.n):r},separators=(',',':'))+'\n')

