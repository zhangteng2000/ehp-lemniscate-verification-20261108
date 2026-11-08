"""Independent replay of all radius bands, normalized source coefficients and collinear moment tests."""
from fractions import Fraction as Q
from pathlib import Path
from math import comb,factorial
from functools import lru_cache
import json
from audit_peeled_line import partition_upper,matrix_variance,gram_check,PLO,PHI,root_floor,upper_sqrt,atan_inv
G=Path(__file__).resolve().parent
scale=1<<140;D=128;nu=Q(1000,433)
def ceil(v):return Q((v*scale).__ceil__(),scale)
def orient(m,j):return comb(m,j) if 0<=j<=m else 0
@lru_cache(None)
def band(m):
    if m==2:return [(Q(1),Q(1),[Q(1),Q(0),Q(1)])]
    previous=[max(c[2][j] for c in band(m-1)) for j in range(m)]
    ts=[Q(0),Q(1,4),Q(1,2),Q(1),Q(2),Q(3),Q(4),Q(5),Q(7),Q(9),Q(12),Q(16),
        Q(m-1,4),Q(m-1,2),Q(3*(m-1),4)]
    shifts=[]
    for ell in range(m):
        v=[Q(comb(m-1-ell,j),(m-1)**j) for j in range(m-ell)]
        shifts.append([abs(a-b) for a,b in zip(v+[0],[0]+v)])
    out=[]
    for cell in range(D):
        a,b=Q(cell,D),Q(cell+1,D);prefs=[Q(0)]
        for j in range(1,m):
            c=min(j*b,nu-a)
            if j+1<m:c=min(c,Q(m+j-1,m-j-1)-a)
            if m>=5:c=min(c,1+j*(nu+1)/(m-1)-a)
            prefs.append(c)
        slopes=partition_upper(prefs)
        exp=[Q(1)]+[Q(0)]*m
        for z in slopes:
            # Absolute exponential coefficients of (-1)^j/j! convolved with the geometric series.
            if z==0:continue
            f=[];s=Q(0)
            for j in range(m+1):
                s+=Q((-1)**j,factorial(j));assert s>=0;f.append(ceil(z**j*s))
            exp=[ceil(sum(exp[j-i]*f[i] for i in range(j+1))) for j in range(m+1)]
        cm=min((1-Q(m,m-1)*a*t/(1+t))/(1-t/(m-1)) for t in ts if t<m-1)
        avg=min((nu-a)/(m-1),b,(nu*cm+b)/(m-1))
        arr=[Q(1),Q(0)];leading=[Q(1),Q(0)]+[Q(j-1,factorial(j))*b**j for j in range(2,m+1)]
        for j in range(2,m+1):
            expv=sum(leading[i]*exp[j-i] for i in range(j+1))
            sh=sum(previous[i]*cm**i*b**(j-i)*shifts[i][j-i] for i in range(min(m-1,j)+1))
            simple=orient(m-1,j)*avg**j+orient(m-1,j-1)*b*avg**(j-1)
            arr.append(ceil(min(expv,sh,simple,comb(m,j)*min(b,nu/m)**j)))
        out.append((a,b,arr))
    return out

@lru_cache(None)
def trig(k):
    lo=PLO/(2*k);hi=PHI/(2*k)
    def t(v,p,d):return sum((-1)**j*v**(2*j+p)/factorial(2*j+p) for j in range(d+1))
    return lo/t(hi,1,12),hi/t(lo,1,11),1-lo*t(hi,0,11)/t(hi,1,12)
def cond(n,k,a,env,T):
    m=n-k;deg=range(2,m+1);q=1/T**2
    v=Q(1) if a==0 else min(Q(1),nu/(m*a))
    # Enclose whole integer powers before taking their square roots.
    bbot=root_floor((1-q*q)**m,2,1<<110);btop=upper_sqrt((1+v)**m)
    atop=min(1/(1-q*q),1+v*q*q/(1-q))
    if a:
        bottom=1-sum(env[j]*(q/a)**j for j in deg)
        bbot=max(bottom,bbot)
        btop=min(btop,1+sum(env[j]/a**j for j in deg))
        atop=min(atop,1+sum(j*env[j]*(q/a)**j for j in deg)/(m*bbot))
    invroot=(2*atop)**k*btop/(bbot*T**(2*n));assert invroot<1
    blo,bhi,w=trig(k)
    cost=2*m*T*(bhi-1+upper_sqrt(atop/2)*w*invroot/(2*n-1))
    return 2*Q(913,1000)*blo,cost
def check(n,L,root,params):
    c=[None]*(n-1)+[Q(0),Q(1)];wx=root['energy']
    assert len(wx)==n-3
    for idx,k in enumerate(range(n-2,1,-1)):
        m=n-k;bb=band(m);ws=list(map(Q,wx[idx]));ww=params[str(n)][str(k)]
        assert len(ws)==len(ww)==len(bb)
        new=[Q(0)]*(m+1)
        for A,ct,(a,b,env) in zip(ws,ww,bb):
            K,E=cond(n,k,a,env,Q(ct,100))
            assert A>0
            if A*(K-b*E)<L:
                assert a>0
                radius=a*A
                total=Q(0)
                for j in range(2,m+1):
                    z=env[j]*A**j
                    if c[n-j] is not None:z=min(z,c[n-j]*Q(comb(m,j),comb(n,j)))
                    total+=z*radius**(m-j)
                assert total<radius**m,(n,k,ct,float(total/radius**m))
            for j in range(2,m+1):
                new[j]=max(new[j],ceil(env[j]*A**j)*Q(comb(n,j),comb(m,j)))
        for j in range(2,m+1):c[n-j]=new[j] if c[n-j] is None else min(new[j],c[n-j])
    cc=root['coefficients']
    assert all(v<=Q(cc[j]) for j,v in enumerate(c) if v is not None)
    return c
def main():
    l,h=atan_inv(5);ll,hh=atan_inv(239);assert PLO<16*l-4*hh<16*h-4*ll<PHI
    pars=json.loads((G/'threshold9_contours.json').read_text());out={}
    vec=json.loads((G/'conditional_line_primary.json').read_text());small=json.loads((G/'conditional_nonic_primary.json').read_text())
    for n in range(9,15):
        row=(vec if n>9 else small)[str(n)]
        c=check(n,Q(7,8) if n>9 else Q(41,50),row,pars)
        if n==9:
            tar=[Q(i,1000) for i in (103,14,84,311,521,785)]
            assert all(c[j]<=tar[j-2] for j in range(2,8))
            c[2:8]=tar;var=matrix_variance(n,c,Q(12,25))/2;assert var<Q(94,100)
            out[n]=[float(var)]
        else:
            u,r=map(Q,row['test'][:2]);sm=matrix_variance(n,c,u)/2;gr=gram_check(n,c,u,r)
            assert max(sm,gr)<Q(96,100);out[n]=[float(sm),float(gr)]
        print(n,out[n],flush=True)
    (G/'conditional_independent.json').write_text(json.dumps(out,indent=1)+'\n')
if __name__=='__main__':assert __debug__;main()

