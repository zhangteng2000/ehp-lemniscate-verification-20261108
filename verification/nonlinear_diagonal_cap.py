"""Upper bounds on Re h_n in the cancelling direction using the exterior area theorem."""
from fractions import Fraction as Q
from math import factorial,comb,isqrt
from collections import Counter
from functools import lru_cache
import argparse
def split(n,k,low=2):
    if k==0:
        if n==0:yield ()
        return
    for j in range(low,n//k+1):
        for tail in split(n-j,k-1,j):yield (j,)+tail
def terms(n):
    M=n-2;res=[]
    for k in range(2,n//2+1):
        d=Q(1,2*k)-Q(comb(2*k,k),4**k*(2*k-1))
        for pp in split(n,k):
            exps=Counter(j-1 for j in pp);v=1;f=1
            for i,q in exps.items():v*=i**q;f*=factorial(q)
            S=1<<40;root=Q(isqrt(v*S*S)+1,S);w=d*factorial(k)*root/f
            res.append((w,exps))
    return res
def verify(n,goal,den=256):
    # Rational box proof on the positive Euclidean unit ball. Last two variables
    # need no symmetry. The expression is monotone separately in each variable.
    t=terms(n);m=max(max(s) for _,s in t);stack=[((0,)*m,(den,)*m,0)]
    counts=[0,0];depth=0
    while stack:
        lo,hi,d=stack.pop();depth=max(d,depth);budget=den*den-sum(x*x for x in lo)
        if budget<0:counts[0]+=1;continue
        hi=tuple(min(v,isqrt(budget+a*a)+1) for a,v in zip(lo,hi))
        # Cap at rational upper endpoints; every feasible point remains inside.
        def term(box,w,exp):
            for i,p in exp.items():w*=Q(box[i-1],den)**p
            return w
        ceiling=sum(term(hi,w,exp) for w,exp in t)
        if ceiling<goal:counts[1]+=1;continue
        ix=max(range(m),key=lambda i:hi[i]-lo[i]);mid=(hi[ix]+lo[ix])//2
        assert lo[ix]<mid<hi[ix],('uncovered',n,lo,hi,float(ceiling))
        A=lo[:ix]+(mid,)+lo[ix+1:];B=hi[:ix]+(mid,)+hi[ix+1:]
        stack.extend([(lo,B,d+1),(A,hi,d+1)])
    return counts,depth
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=9);p.add_argument('--bound',type=Q,default=Q(9,10))
    a=p.parse_args();print(a.n,a.bound,verify(a.n,a.bound))

