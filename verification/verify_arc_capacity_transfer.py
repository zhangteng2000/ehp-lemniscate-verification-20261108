"""Exact all-phase circle reductions in degrees 6 and 7 using the capacity of the marked arc."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
from math import factorial
from functools import lru_cache
from audit_peeled_line import root_floor,upper_sqrt,PLO,PHI
from audit_low_degree_transfer import bino,const
from audit_quintic_rational import I,taylor
from diagonal_curvature_bound import certified
from coherent_projection_transfer import constants
from nonlinear_diagonal_cap import verify
import json,argparse
G=Path(__file__).resolve().parent;S=1<<48
def up(y):return Q((S*y).__ceil__(),S)
def dn(y):return -up(-y)
@lru_cache(None)
def raw(n,M):
    x,b,il,ih,gap,H=const(n);rows=[]
    for m in range(1,M+1):
        P=[];cur=Q(0)
        for i in range(m):cur+=(-1)**m*bino(b,i)*bino(-b,m-i);P.append(cur)
        ratios=[Q(1)]
        for k in range(1,m):ratios.append(ratios[-1]*Q(k-1-b,k+b))
        A=[sum(P[i]*ratios[abs(k-i)] for i in range(m)) for k in range(m)]
        def ip(v):
            z=sum(vv*ww for vv,ww in zip(A,v));p=P[0]*v[0]
            return I(dn(min(il*z-p,ih*z-p)),up(max(il*z-p,ih*z-p)))
        norm=ip(P);assert norm.l>0;group={};ones=[]
        from math import comb
        for j in range(1,m+1):
            for l in range(m+1-j):
                v=[Q(0)]*m;r=m-j-l
                for i in range(r+1):
                    h=(-1)**r*bino(b-j-l,i)*bino(-b-j,r-i)
                    for e in range(j):v[i+l+e]+=h*(-1)**e*comb(j-1,e)
                a=ip(v)
                if j==1:ones.append(a)
                else:
                    cap=H*Q((j+l)**(j+l),j**j*(l**l if l else 1))
                    group[j+l]=up(group.get(j+l,0)+a.absupper()*cap)
        rows.append((norm,ones,group))
    return rows
def response(x,c0,c1,H=24):
    # Expand at the exact midpoint with a uniform derivative bound on [-1,1].
    c=(c0+c1)/2;w=(c1-c0)/2;lam=x/2
    a=Q(1);b=x*c;total=Q(1)+max(Q(0),abs(b)-x*w)**2;bound=x
    for k in range(1,H):
        a,b=b,(2*(k+lam)*c*b-(k-1+2*lam)*a)/(k+1)
        bound*=Q(k+1+x,k)
        total+=max(Q(0),abs(b)-bound*w)**2
    return dn(total)
def run(n,L,M,t_override=None,signed=False):
    cap={3:Q(0),6:Q(27,100),7:Q(36,100)}[n]
    cc=constants(n);x=Q(1,n);t={3:Q(1,110),6:Q(1,230),7:Q(1,392)}[n]
    if t_override is not None:t=Q(t_override)
    hn=cap/2+Q(n-1,2*n)*t if signed else None
    small=certified(n,3,t,Q(2,5),cc,hn_bound=hn);assert small<L
    if n>3:
        from audit_capacity_majorants_exact import verify as verify_majorants
        verify_majorants()
    x,b,il,ih,gap,H=const(n);rows=raw(n,M)
    tt=2*t/(1+t*t);end=Q(173,200);assert end**2<Q(3,4)
    todo=[(tt,end),(end,Q(1))];leaves=[];rho=[Q(k,40) for k in range(1,40)]
    while todo:
        l,h=todo.pop();c0=dn(root_floor(1-h*h,2));c1=min(Q(1),up(upper_sqrt(1-l*l)))
        jmin=max(cc[-1].l,response(x,c0,c1));gain=Q(0)
        if h<=end:
            ac=up(1/root_floor(2*(1+c0),2));rr=up(h/(1+c0))
            re=[]
            for v in range(M):
                er=cap*(ac if v==0 else min(min(1/(1-r*r),ac/(1-r)**2)/r**v for r in rho))
                if v%2:
                    p0=Q(2*(-1)**((v+1)//2),2*n);d=Q(n-3,2*n)
                    lo,hi=p0-d,p0+d
                elif v:
                    hi=Q(n-3,2*n)*min(1,(v+1)*rr);lo=-hi
                else:lo=0;hi=Q(n-3,2*n)*rr
                re.append(I(dn(lo-er),up(hi+er)))
            ch=[Q(1)]+[up(c1**j) for j in range(1,M+1)];hp=[Q(1)]+[up(h**j) for j in range(1,M+1)]
            for m,(nm,ones,grp) in enumerate(rows,1):
                v=dn(c0**m*nm.l)
                for r,z in grp.items():v-=up(z*hp[r]*ch[m-r])
                for k,a in enumerate(ones):
                    z=(a*re[k]).l
                    if z<0:v+=dn(z*hp[k+1]*ch[m-k-1])
                    else:v+=dn(z*l**(k+1)*c0**(m-k-1))
                gain+=dn(max(0,v)**2/nm.h)
        rem=up((gap-2*PLO*gain)/(L*jmin))
        if rem<=0 or rem**n<l:leaves.append((l,h));continue
        if h-l<Q(1,1<<25):raise ArithmeticError((n,float(l),float(h),float(rem)/float(l)**float(x)*float(L)))
        c=(l+h)/2;todo.extend(((l,c),(c,h)))
    leaves.sort();assert leaves[0][0]==tt and leaves[-1][1]==1 and all(i[1]==j[0] for i,j in zip(leaves,leaves[1:]))
    ans={'degree':n,'conditional_line':str(L),'nonlinear_cap':str(cap),'rows':M,'tau':str(t),
         'signed_first_term':signed,'curvature_upper_display':float(small),'intervals':[[str(u),str(v)] for u,v in leaves]}
    (G/f'capacity_transfer_{n}.json').write_text(json.dumps(ans,indent=1)+'\n')
    print(n,L,M,len(leaves),'closed diameter intervals')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--L',type=Q,required=True)
    p.add_argument('--M',type=int,required=True)
    p.add_argument('--tau',type=Q);p.add_argument('--signed',action='store_true')
    a=p.parse_args();assert __debug__;run(a.n,a.L,a.M,a.tau,a.signed)

