"""Exact sufficient cubic areas for a full sextic or septic EHP proof. No cubic cover is assumed here."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import argparse,json
from seeded_coefficient_descent import descend
from audit_peeled_line import matrix_variance,PLO,PHI
G=Path(__file__).resolve().parent;S=1<<38
def floor(a):return (a*S).__floor__()
def ceil(a):return (a*S).__ceil__()
@lru_cache(None)
def sine(a,b):
    d=Q(0) if a.__ceil__()<=b else min(a-a.__floor__(),b.__ceil__()-b)
    z=PLO*d/2;w=PHI*d/2
    return floor(2*(z-w**3/6)*(1-w*w/2+z**4/24-w**6/720))
@lru_cache(None)
def vals(n,nt,ng):
    return [[(sine(Q(n*i,nt),Q(n*(i+1),nt)),
              sine(Q(i,nt)+Q(g,2*n*ng),Q(i+1,nt)+Q(g+1,2*n*ng))) for i in range(nt)] for g in range(ng)]
def strips(n,c,u,nt=432,ng=8,nh=512):
    m=n-1;a=vals(n,nt,ng);b=[Q(0)]+[Q(i+1,n)*c[i+1] for i in range(1,m+1)]
    coef=[]
    for k in range(nh+1):
        h=Q(k,nh);v=w=z=Q(0)
        for j in range(m,-1,-1):
            z=z*h+2*w;w=w*h+v;v=v*h+b[j]
        coef.append((ceil(v),ceil(w),ceil(z)))
    M=m**m
    T={}
    for x in set(x for r in a for x,y in r):
        # T/S^2 majorizes m H^2(m-(x/S)^2)+2 B_(m-2)*(x/S).
        p=Q(m*(m*S*S-x*x),nh*nh);q=2*b[m-2]*x*S
        T[x]=[((p*k*k+q).__ceil__()) for k in range(nh+1)]
    @lru_cache(None)
    def height(x,y):
        if x==0 or y==0:return 0
        E=(u*y/n).__floor__();xx=x*x;e=M*xx;f=M*x**(2*m-2);h=m*m*xx
        def yes(k):
            v,w,z=coef[k];d=E-v;j=T[x][k]
            if d<=0:return False
            r=w*w+z*d
            return e*d**(2*m-2)>r**m or f*d*d>j**m or h*d*d>j*r
        l,r=0,nh
        while r-l>1:
            k=(l+r)//2
            if yes(k):l=k
            else:r=k
        assert l==0 or yes(l)
        return l
    return min(2*PLO*sum(height(x,y) for x,y in ar)/(nt*nh) for ar in a)
def row(n,L,chart,z,q,params):
    c=[None]*(n-1)+[Q(0),Q(1)]
    if chart=='A':f=(L/q)**2;g=z*(L/q)**3
    else:f=z*(L/q)**2;g=(L/q)**3
    c[n-2]=Q(comb(n,2),3)*f;c[n-3]=comb(n,3)*g
    for _ in range(2):descend(n,L,params,c)
    # Bound source moduli by rational multiples of 1/10000. Round upward.
    for i in range(2,n):c[i]=Q((c[i]*10000).__ceil__(),10000)
    l,h=0,2048
    while h-l>1:
        v=(l+h)//2
        if matrix_variance(n,c,Q(v,512))<2:l=v
        else:h=v
    u=Q(l,512)
    if u==0:return False,None
    vv=matrix_variance(n,c,u)/2;assert vv<1
    area=strips(n,c,u)
    return area>L,[str(u),[str(t) for t in c[2:]],str(vv),str(area)]
def generate(n,L):
    records=json.loads((G/f'cubic_shape_requested_{n}.json').read_text())
    params=json.loads((G/'conditional_contours.json').read_text());out=[]
    for a,i,_,v in records:
        z=Q(i+1,32)*(Q(5,16) if a=='A' else Q(9,4));q=Q((100*v).__ceil__(),100)
        for _ in range(15):
            ok,cert=row(n,L,a,z,q,params)
            if ok:break
            q+=Q(1,100)
        assert ok,(n,a,i,q)
        out.append([a,i,str(z),str(q),*cert])
        print(n,a,i,str(q),float(Q(cert[-1])),flush=True) if i%8==7 else None
    assert [(d[0],d[1]) for d in out]==[(k,j) for k in 'AB' for j in range(32)]
    (G/f'proved_cubic_sufficient_profile_{n}.json').write_text(json.dumps(out,indent=1)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--L',type=Q,required=True)
    a=p.parse_args();assert __debug__;generate(a.n,a.L)

