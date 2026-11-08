"""Exact derivative-radius and collinear-line certificate. Run with standard Python -S -B."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial,isqrt
from pathlib import Path
import json,argparse
from verify_peeled_line import choose,slopes,variance,gram
from verify_interpolation_degrees import LOWPI,HIGHPI,pi_check
G=Path(__file__).resolve().parent
N=128; B=1<<58; C=Q(913,1000); NU=Q(1000,433)
def up(x):
    x=Q(x);return Q(-((-x.numerator*B)//x.denominator),B)
def down(x):return -up(-x)
def usqrt(x):
    x=Q(x);i=isqrt(x.numerator*B*B//x.denominator)
    return Q(i+int(i*i*x.denominator!=x.numerator*B*B),B)
def lsqrt(x):return max(Q(0),usqrt(x)-Q(1,B))
def hp(x,m,upper):
    a=x**(m//2)
    if m%2:a*=usqrt(x) if upper else lsqrt(x)
    return up(a) if upper else down(a)
def make_envelope(m,parts=N):
    if m==2:return ((Q(1),Q(1),(Q(1),Q(0),Q(1))),)
    old=universal(m-1,parts);w=[Q(i) for i in (0,1,2,3,4,5,7,9,12,16)]
    w.extend((Q(1,4),Q(1,2),Q(m-1,4),Q(m-1,2),Q(3*(m-1),4)))
    rows=[]
    for it in range(parts):
        l,h=Q(it,parts),Q(it+1,parts);pref=[Q(0)]
        for i in range(1,m):
            b=min(i*h,NU-l)
            if i<m-1:b=min(b,1+Q(2*i,m-i-1)-l)
            if m>=5:b=min(b,1+i*(NU+1)/(m-1)-l)
            pref.append(b)
        lam=slopes(pref);s=[0]+[up(sum(x**j for x in lam)) for j in range(1,m+1)]
        H=[Q(1),Q(0)];d=[Q(1),Q(0)]
        z=min((1-Q(m,m-1)*l*t/(1+t))/(1-t/(m-1)) for t in w if t<m-1)
        assert z>=0
        v=min(h,(NU-l)/(m-1),(NU*z+h)/(m-1));ee=[Q(1),Q(0)]
        for j in range(2,m+1):
            H.append(up(sum(s[t]*H[j-t] for t in range(2,j+1))/j))
            d.append(Q(j-1,factorial(j))*h**j)
            t1=sum(H[j-i]*d[i] for i in range(j+1))
            t2=sum(old[i]*z**i*(h/(m-1))**(j-i)*abs(choose(m-1-i,j-i)
                 -(m-1)*choose(m-1-i,j-i-1)) for i in range(min(m-1,j)+1))
            t3=choose(m-1,j)*v**j+h*choose(m-1,j-1)*v**(j-1)
            ee.append(up(min(t1,t2,t3,comb(m,j)*min(NU/m,h)**j)))
        rows.append((l,h,tuple(ee)))
    return tuple(rows)
make_envelope=lru_cache(None)(make_envelope)
@lru_cache(None)
def universal(m,parts):
    return [max(r[2][j] for r in make_envelope(m,parts)) for j in range(m+1)]
@lru_cache(None)
def trig(k):
    x=LOWPI/(2*k);y=HIGHPI/(2*k)
    def pol(t,par,end):return sum((-1)**s*t**(2*s+par)/factorial(2*s+par) for s in range(end+1))
    low=pol(x,1,7);high=pol(y,1,8)
    return down(x/high),up(y/low),up(1-x*pol(y,0,7)/high)

def errors(n,k,row,tc):
    m=n-k;l,h,e=row;T=Q(tc,100);p=1/T**2;assert T>1
    sig=Q(1) if not l else min(Q(1),NU/(m*l))
    z=hp(1-p*p,m,False);w=hp(1+sig,m,True)
    a=min(1/(1-p*p),1+sig*p*p/(1-p))
    if l:
        qs=[up(e[j]*(p/l)**j) for j in range(2,m+1)]
        z=max(z,1-sum(qs))
        w=min(w,1+sum(up(e[j]/l**j) for j in range(2,m+1)))
        a=min(a,1+sum(j*q for j,q in enumerate(qs,2))/(m*z))
    a=up(a);rho=up((2*a)**k*w/(T**(2*n)*z))
    be_low,be_high,bk=trig(k)
    if rho>=1:return None
    E=up(2*m*T*(be_high-1+usqrt(a/2)*bk*rho/(2*n-1)))
    return down(2*C*be_low),E

def get_caps(n,L,params,parts=N):
    co=[None]*(n-1)+[Q(0),Q(1)];records=[]
    for k in range(n-2,1,-1):
        m=n-k;env=make_envelope(m,parts)
        fac=[Q(comb(n,j),comb(m,j)) for j in range(m+1)]
        q=[None if co[n-j] is None else co[n-j]/fac[j] for j in range(m+1)]
        new=[Q(0)]*(m+1);ws=[]
        ts=params[str(n)][str(k)];assert len(ts)==len(env)
        for (l,h,e),tc in zip(env,ts):
            ke=errors(n,k,(l,h,e),tc);assert ke is not None,(n,k,tc)
            K,E=ke
            geo=up(L/(K-E*h)) if K>E*h else None
            def fits(A):
                if not l:return False
                return sum((e[j] if q[j] is None else min(e[j],q[j]/A**j))/l**j for j in range(2,m+1))<1
            A=geo
            hi=min(geo,Q(64)) if geo is not None else Q(64)
            if fits(hi):
                lo=Q(0)
                for _ in range(16):
                    b=(lo+hi)/2
                    if fits(b):hi=b
                    else:lo=b
                A=up(hi)
            assert A is not None and A>0 and (A*(K-E*h)>=L or fits(A)),(n,k,l)
            ws.append(str(A))
            for j in range(2,m+1):
                y=up(e[j]*A**j*fac[j]);new[j]=max(new[j],y)
        for j in range(2,m+1):
            v=new[j];co[n-j]=v if co[n-j] is None else min(v,co[n-j])
        records.append(ws)
    return co,records

PARAM={10:(Q(66,100),Q(1445,1000)),11:(Q(58,100),Q(1385,1000)),12:(Q(52,100),Q(134,100)),
       13:(Q(45,100),Q(129,100)),14:(Q(39,100),Q(126,100))}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=list(PARAM))
    p.add_argument('--L',type=Q,default=Q(7,8));p.add_argument('--parts',type=int,default=N)
    p.add_argument('--parameters',default='threshold9_contours.json');p.add_argument('--out',default='conditional_line_primary.json')
    opt=p.parse_args();assert __debug__;pi_check()
    pars=json.loads((G/opt.parameters).read_text());out={}
    for n in opt.degrees:
        co,ws=get_caps(n,opt.L,pars,opt.parts)
        r={'coefficients':[None if v is None else str(v) for v in co],'energy':ws}
        if n in PARAM and opt.L==Q(7,8):
            u,rad=PARAM[n];sm=variance(n,co,u)/2;lg=gram(n,co,u,rad)
            assert sm<Q(96,100) and lg<Q(96,100),(n,float(sm),float(lg))
            r['test']=[str(u),str(rad),float(sm),float(lg)]
            print(n,'variance/2 and Gram',float(sm),float(lg),flush=True)
        else:print(n,[None if v is None else round(float(v),7) for v in co],flush=True)
        out[n]=r
    (G/opt.out).write_text(json.dumps(out,separators=(',',':'))+'\n')

