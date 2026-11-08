"""Rational circle transfer with the signed first exterior row in degrees six, seven and eight."""
from fractions import Fraction as Q
from math import comb
from functools import lru_cache
from pathlib import Path
import argparse,json
from coherent_projection_transfer import constants,rows
from diagonal_curvature_bound import certified
from audit_quintic_rational import I,PI
from nonlinear_diagonal_cap import verify

# Where there is more than one row use a positive lower bound on each row individually.
degree={6:(Q(19,25),Q(7,25),Q(1,230),18),7:(Q(19,25),Q(37,100),Q(1,392),10),
        8:(Q(39,50),Q(53,100),Q(1,512),7)}
def certify(n,params=None,use_jmin=True):
    L,K,t,M=degree[n] if params is None else params
    const=constants(n);x,b,C,I0,U,J=const;gap=(U-C).h
    cap=verify(n,K)
    curve=certified(n,3,t,Q(2,5),const)
    assert curve<L,('short curvature',n,float(curve),float(L))
    rows0=rows(n,M,const);start=2*t/(1+t*t)
    stack=[(start,Q(3,4)),(Q(3,4),Q(1))];done=[]
    while stack:
        lo,hi=stack.pop();extra=Q(0)
        if hi<=Q(3,4):
            c=(I(1-hi*hi)**Q(1,2)).l;cb=(I(1-lo*lo)**Q(1,2)).h
            for m,(N,band) in enumerate(rows0,1):
                err=sum(a*hi**j*cb**(m-j) for j,a in band.items())
                v=max(Q(0),c**m*N.l-err)**2/N.h
                if m==1:v=max(v,(I0.l-1)*max(Q(0),b*c-K*hi)**2)
                extra+=v
        target=(gap-2*PI.l*extra)/(L*(J.l if use_jmin else Q(1)))
        if target<=0 or target**n<lo:done.append([lo,hi]);continue
        assert hi-lo>Q(1,2**30),('diameter',n,float(lo),float(hi))
        mid=(lo+hi)/2;stack.extend([(lo,mid),(mid,hi)])
    done.sort();assert done[0][0]==start and done[-1][1]==1
    assert all(j[1]==k[0] for j,k in zip(done,done[1:]))
    return {'L':str(L),'K':str(K),'tau0':str(t),'M':M,'j_min':use_jmin,'curvature_display':float(curve),
            'ball_check':cap,'intervals':[[str(a),str(b)] for a,b in done]}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[6,7,8]);p.add_argument('--pure-jensen',action='store_true');a=p.parse_args()
    for n in a.n:
        r=certify(n,use_jmin=not a.pure_jensen);print(n,len(r['intervals']),r['curvature_display'],flush=True)
        Path(__file__).with_name(f'new_signed_circle_{n}.json').write_text(json.dumps(r,indent=1)+'\n')

