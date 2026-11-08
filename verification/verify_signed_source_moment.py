"""Certified signed top moment of a real slice, with full coverage of both continuous phases."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from math import isqrt,comb,factorial
from pathlib import Path
from bisect import bisect_right
import argparse,json
from audit_peeled_line import PLO,PHI,matrix_variance
from seeded_coefficient_descent import descend
G=Path(__file__).resolve().parent;S=1<<30
def up(x):return x.__ceil__()
@lru_cache(None)
def sine_at(t):
    t=t%2;neg=t>=1;t=t%1;d=min(t,1-t)
    if d==0:return (0,0)
    if d==Q(1,2):return ((-S,-S) if neg else (S,S))
    a=PLO*d;b=PHI*d
    l=sum((a**(2*j+1) if j%2==0 else -b**(2*j+1))/factorial(2*j+1) for j in range(8))
    h=sum((b**(2*j+1) if j%2==0 else -a**(2*j+1))/factorial(2*j+1) for j in range(9))
    lo=(S*l).__floor__();hi=up(S*h)
    return (-hi,-lo) if neg else (lo,hi)
@lru_cache(None)
def sinebox(a,b):
    lo=min(sine_at(a)[0],sine_at(b)[0]);hi=max(sine_at(a)[1],sine_at(b)[1])
    if ((a-Q(1,2))/2).__ceil__()<=((b-Q(1,2))/2).__floor__():hi=S
    if ((a-Q(3,2))/2).__ceil__()<=((b-Q(3,2))/2).__floor__():lo=-S
    return lo,hi
@lru_cache(None)
def angle_sets(n,d,nt,np0,ng):
    data=[]
    for j in range(nt):
        a,b=Q(j,nt),Q(j+1,nt);v,w=sinebox(n*a,n*b);x=max(0,v,-w);ymax=max(w,-v)
        if x==0:continue
        ys=[max(0,*(lambda p:(p[0],-p[1]))(sinebox(d*a+Q(k,ng),d*b+Q(k+1,ng)))) for k in range(ng)]
        signed=[]
        for p in range(np0):
            y,z=sinebox((n-2)*a+Q(p,n*np0),(n-2)*b+Q(p+1,n*np0))
            signed.append(y if v>0 else -z)
        data.append((x,ymax,signed,ys))
    return data
@lru_cache(None)
def base_table(m,x,nh,hi):
    return [-( -m*k*k*(m*S*S-x*x)//(nh*nh)) for k in range(hi+1)]
def coverage(n,d,c,u,toplo,tophi,mesh,paired=False):
    nt,np0,ng,nh,hi=mesh;m=n-d;C=comb(n,d);B=Q(comb(n-2,d),C)
    geom=angle_sets(n,d,nt,np0,ng);V=[];VE=[]
    for k in range(hi+1):
        h=Q(k,nh);v=e=Q(0)
        for j in range(m,0,-1):
            b=Q(comb(j+d,d),C)*c[j+d];v=h*v+b;e=h*e+(b if j%2==0 else 0)
        V.append(up(S*h*v));VE.append(up(S*h*e))
    T0=[-(-m*(m-1)*k*k*S*S//(nh*nh)) for k in range(hi+1)]
    acc=[[0]*ng for p in range(np0)];unit=m**m
    for x,xhi,top,ys in geom:
        power=unit*x**(2*m-2);tbase=base_table(m,x,nh,hi)
        for p,b in enumerate(top):
            offset=up(-2*B*b*(x*toplo if b>=0 else xhi*tophi));thresholds=[];one=[]
            for k in range(1,hi+1):
                t=tbase[k]+offset
                if T0[k]+offset<0:thresholds.append(0);one.append(0);continue
                if u==0:break
                assert t>=0
                root=isqrt(t**m//power)+1;need=V[k]+root
                least=-( -C*need*u.denominator//u.numerator)
                if paired:
                    least_one=-(-C*(VE[k]+root)*u.denominator//u.numerator)
                    if least_one>S:break
                    one.append(least_one)
                elif least>S:break
                thresholds.append(least)
            assert all(a<=b for a,b in zip(thresholds,thresholds[1:]))
            assert all(a<=b for a,b in zip(one,one[1:]))
            sums=acc[p]
            for j,y in enumerate(ys):sums[j]+=bisect_right(thresholds,y)+(bisect_right(one,y) if paired else 0)
    num=min(map(min,acc))
    return (1 if paired else 2)*PLO*num/(nt*nh)
def test(n,L,kind,i,q,mesh,paired=False,parts=2,conditioned=False):
    P=32;ztop=Q(5,16) if kind=='A' else Q(9,4)
    z0=i*ztop/P;z1=(i+1)*ztop/P;Smax=L/q
    pars=json.loads((G/'conditional_contours.json').read_text());rec=[];closed=True
    cuts=[Q(0),Q(4,5),Q(9,10),Q(19,20),Q(1)]
    if conditioned:cuts=[Q(j,100) for j in (0,80,88,92,95,97,99,100)]
    for l,h in zip(cuts,cuts[1:]):
        sl=l*Smax;su=h*Smax
        flo=sl**2*(1 if kind=='A' else z0);f=su**2*(1 if kind=='A' else z1)
        g=su**3*(z1 if kind=='A' else 1);toplo=Q(comb(n,2),3)*flo
        c=[None]*(n-1)+[Q(0),Q(1)];c[n-2]=Q(comb(n,2),3)*f;c[n-3]=comb(n,3)*g
        lower=[Q(0)]*(n+1);lower[n-2]=toplo;lower[n-3]=comb(n,3)*sl**3*(z0 if kind=='A' else 1)
        good=True
        for _ in range(2):
            if descend(n,L,pars,c,lower=lower if conditioned else None) is None:good=False;break
        if not good:continue
        for j in range(2,n-1):c[j]=Q(up(10000*c[j]),10000)
        if toplo>c[n-2]:continue
        for part in range(parts):
            cc=c.copy();v=part*cc[2]/parts;cc[2]=(part+1)*cc[2]/parts
            if conditioned:
                bounded=lower.copy();bounded[2]=v
                if descend(n,L,pars,cc,lower=bounded) is None:continue
                if any(cc[j]<bounded[j] for j in range(2,n-1)):continue
            a,b=0,4096
            while b-a>1:
                mid=(a+b)//2
                if matrix_variance(n,cc,Q(mid,1024))<2:a=mid
                else:b=mid
            u=Q(a,1024);ans=Q(0)
            if u and matrix_variance(n,cc,u)<2:
                ans=coverage(n,1,cc,u,toplo,c[n-2],mesh,paired)
            if ans<=L:
                # With the below split on c2 this argument is unconditional in c1.
                ans=max(ans,coverage(n,2,cc,v,toplo,c[n-2],mesh,paired))
            rec.append([str(l),str(h),part,str(u),str(ans)])
            if ans<=L:closed=False
    return closed,rec
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--L',type=Q,required=True)
    p.add_argument('--kind',default='B');p.add_argument('--q',type=Q,required=True);p.add_argument('--bins',nargs='+',type=int,required=True)
    p.add_argument('--mesh',nargs=5,type=int,default=[420,8,96,400,200]);p.add_argument('--paired',action='store_true');p.add_argument('--conditioned',action='store_true');p.add_argument('--parts',type=int,default=2);a=p.parse_args();assert __debug__
    for i in a.bins:
        ok,r=test(a.n,a.L,a.kind,i,a.q,tuple(a.mesh),a.paired,a.parts,a.conditioned);shape=[a.n,a.kind,i,str(a.q),a.mesh,ok,r]
        if a.parts!=2:shape.append(a.parts)
        tail=('' if a.parts==2 else f'_p{a.parts}')+('_cond' if a.conditioned else '')
        dest=G/f'{"paired" if a.paired else "signed"}_moment_{a.n}_{a.kind}_{i}_{a.q.numerator}_{a.q.denominator}{tail}.json'
        dest.write_text(json.dumps(shape,indent=1)+'\n')
        print(a.n,a.kind,i,a.q,'closed',ok,'min area',float(min(Q(v[-1]) for v in r)),flush=True)

