"""Independent sine enclosures and integral-height witnesses for the remaining degree reductions."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from bisect import bisect_right
from math import comb,isqrt
from audit_peeled_line import PLO,PHI
G=Path(__file__).resolve().parent;D=1<<34
@lru_cache(None)
def table():
    d=(G/'quintic_sine_independent.txt').read_text().splitlines();assert d[0]=='4096 46' and len(d)==8193
    return [tuple(map(int,x.split())) for x in d[1:]]
@lru_cache(None)
def point(a):
    # Linear interpolation error for sin(pi t) is <= pi^2/(8*4096^2).
    v=a*4096;j=v.__floor__();w=v-j;e=Q(10,8*4096**2) if w else 0;f=table()
    l=((1-w)*f[j%8192][0]+w*f[(j+1)%8192][0])/(1<<46)-e
    h=((1-w)*f[j%8192][1]+w*f[(j+1)%8192][1])/(1<<46)+e
    return (l*D).__floor__(),(h*D).__ceil__()
@lru_cache(None)
def sn(a,b):
    l=min(point(a)[0],point(b)[0]);h=max(point(a)[1],point(b)[1])
    for z,v in ((Q(1,2),1),(Q(3,2),-1)):
        if ((a-z)/2).__ceil__()<=((b-z)/2).__floor__():
            if v==1:h=D
            else:l=-D
    return l,h
@lru_cache(None)
def geo(n,d,K,A,P,typ):
    grid=[]
    for j in range(K):
        l,h=Q(j,K),Q(j+1,K);xx,yy=sn(n*l,n*h);x=max(0,xx,-yy)
        if x==0:continue
        phases=[]
        for k in range(P):
            span=Q(1,2*n) if typ=='unsigned' else Q(1)
            zl,zh=sn(d*l+k*span/P,d*h+(k+1)*span/P);phases.append(max(0,zl,-zh))
        if typ=='unsigned':sg=[]
        else:
            sg=[]
            for i in range(A):
                aa,bb=sn((n-2)*l+Q(i,A*n),(n-2)*h+Q(i+1,A*n))
                sg.append(aa if xx>0 else -bb)
        grid.append((x,max(yy,-xx),phases,sg))
    return grid
def source(n,d,c,u,lo,hi,typ,mes):
    K,A,P,N,Hmax=mes;C=comb(n,d);m=n-d;rows=geo(n,d,K,A,P,'unsigned' if typ=='unsigned' else 'signed')
    # Exact rational Horner evaluations, followed by separate upward dyadic enclosures.
    vals=[]
    for k in range(Hmax+1):
        r=Q(k,N);v=e=w=z=Q(0)
        for j in range(m,-1,-1):
            coeff=Q(comb(d+j,d),C)*c[d+j] if j else Q(0)
            z=z*r+2*w;w=w*r+v;v=v*r+coeff;e=e*r+(coeff if j%2==0 else 0)
        vals.append(tuple((a*D).__ceil__() for a in (v,e,w,z)))
    if typ=='unsigned':
        def get(x,y):
            cached=(x,y)
            if cached in heights:return heights[cached]
            cv=(u*y/C).__floor__();h0,h1=0,Hmax+1;f=m**m;fac=f*x**(2*m-2)
            def yes(k):
                v,_,w,z=vals[k];b=cv-v
                if b<=0:return False
                T=(Q(m*k*k*(m*D*D-x*x),N*N)+2*Q(comb(n-2,d),C)*hi*x*D).__ceil__()
                V=w*w+b*z
                return m**m*x*x*b**(2*m-2)>V**m or fac*b*b>T**m or m*m*x*x*b*b>T*V
            while h1-h0>1:
                k=(h1+h0)//2
                if yes(k):h0=k
                else:h1=k
            assert h0==0 or yes(h0);heights[cached]=h0;return h0
        heights={};acc=[0]*P
        for x,xu,ys,_ in rows:
            for j,y in enumerate(ys):acc[j]+=get(x,y)
        return 2*PLO*min(acc)/(K*N)
    # Here integer thresholds for uY/C are recovered from the real-root product
    # with exact division; the argument of the sine supplies the signed top.
    acc=[[0]*P for _ in range(A)];signedQ=2*Q(comb(n-2,d),C)
    neg=[Q(m*(m-1)*k*k*D**2,N**2).__ceil__() for k in range(Hmax+1)]
    bases={}
    for x,xu,ys,sg in rows:
        if x not in bases:bases[x]=[Q(m*k*k*(m*D*D-x*x),N**2).__ceil__() for k in range(Hmax+1)]
        base=bases[x];den=m**m*x**(2*m-2)
        for a,b in enumerate(sg):
            off=(-signedQ*(hi*xu if b<0 else lo*x)*b).__ceil__()
            test=[];odd=[]
            for k in range(1,Hmax+1):
                if neg[k]+off<0:test.append(0);odd.append(0);continue
                if u<=0:break
                upper=base[k]+off;assert upper>=0
                root=isqrt(upper**m//den)+1
                V,E,_,_=vals[k]
                loY=-(-C*(root+(E if typ=='paired' else V))*u.denominator//u.numerator)
                if loY>D:break
                odd.append(loY);test.append(-(-C*(root+V)*u.denominator//u.numerator))
            # Prefix sums over every accepted angular phase at the corresponding height.
            arr=acc[a]
            for p,v in enumerate(ys):arr[p]+=bisect_right(test,v)+(bisect_right(odd,v) if typ=='paired' else 0)
    return (1 if typ=='paired' else 2)*PLO*min(map(min,acc))/(K*N)

