"""Verify universal coefficient caps, scalar witnesses and a companion-matrix variance independently."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb,factorial,isqrt
from pathlib import Path
import argparse,json
G=Path(__file__).resolve().parent
UNIT=1<<250;PLO=Q("3.1415926535897932384626433832795028");PHI=PLO+Q(1,10**34)
L=Q(49,60);C0=Q(913,1000);NU=Q(1000,433)
def up(t):
    t=Q(t);return Q(-((-t.numerator*UNIT)//t.denominator),UNIT)
def root_floor(t,n,scale=1<<76):
    t=Q(t);ln=0;hn=scale;v=t.numerator*scale**n
    while hn**n*t.denominator<=v:hn*=2
    while hn-ln>1:
        mid=(ln+hn)//2
        if mid**n*t.denominator<=v:ln=mid
        else:hn=mid
    return Q(ln,scale)
def upper_sqrt(t):return root_floor(t,2)+Q(1,1<<76)
def atan_inv(d):
    a=sum(Q((-1)**j,(2*j+1)*d**(2*j+1)) for j in range(50))
    return a,a+Q(1,101*d**101)
def e_bound(n,k):
    x=PHI/(2*k);lo=PLO/(2*k)
    def trig(v,par,N):
        return sum((-1)**j*v**(2*j+par)/factorial(2*j+par) for j in range(N+1))
    A=1/(k*trig(lo,1,9))-2/PHI
    B=2/PLO-trig(x,0,9)/(k*trig(x,1,10))
    return up(2*PHI*(n-k)*(A+4*B*Q(95,117)**n/(2*n-1)))

def partition_upper(prefix):
    # Select the highest future slope at each new anchor and verify the resulting majorant.
    i=0;ans=[]
    while i<len(prefix)-1:
        h=max(range(i+1,len(prefix)),key=lambda j:((prefix[j]-prefix[i])/(j-i),j))
        ans.extend([(prefix[h]-prefix[i])/(h-i)]*(h-i));i=h
    assert all(a>=b>=0 for a,b in zip(ans,ans[1:]))
    assert all(sum(ans[:i])>=p for i,p in enumerate(prefix))
    return ans

@lru_cache(None)
def caps_by_degree(m):
    if m==2:return [Q(1),Q(0),Q(1)]
    prior=caps_by_degree(m-1)
    A=[Q(1),Q(0)]+[Q(0)]*(m-1);D=(1+NU)/(m-1)
    ts=[Q(i) for i in (0,1,2,3,4,5,7,9,12,16)]+[Q(1,4),Q(1,2),Q(m-1,4),Q(m-1,2),Q(3*(m-1),4)]
    shifted=[]
    for ell in range(m):
        degree=m-1-ell
        base=[Q(comb(degree,j),(m-1)**j) for j in range(degree+1)]
        shifted.append([abs(a-b) for a,b in zip(base+[0],[0]+base)])
    for slot in range(24):
        u=Q(slot,24);v=Q(slot+1,24)
        pref=[Q(0)]
        for i in range(1,m):
            constraints=[i*v,NU-u]
            if i<m-1:constraints.append(Q(m+i-1,m-i-1)-u)
            if m>=5:constraints.append(1+i*D-u)
            pref.append(min(constraints))
        lam=partition_upper(pref)
        # Multiply the positive factors exp(-lambda x)/(1-lambda x).
        h=[Q(1)]+[Q(0)]*m
        for val in lam:
            if not val:continue
            ex=[Q(1)];part=Q(1);cur=Q(1)
            for j in range(1,m+1):
                cur*=-Q(1,j);part+=cur
                assert part>=0
                ex.append(up(part*val**j))
            h=[up(sum(ex[q]*h[j-q] for q in range(j+1))) for j in range(m+1)]
        exact_root=[Q(1),Q(0)]+[Q(j-1,factorial(j))*v**j for j in range(2,m+1)]
        AA=min((1-Q(m,m-1)*u*t/(1+t))/(1-t/(m-1)) for t in ts if t<m-1)
        for j in range(2,m+1):
            one=sum(h[q]*exact_root[j-q] for q in range(j+1))
            two=sum(prior[ell]*AA**ell*shifted[ell][j-ell]*v**(j-ell)
                    for ell in range(min(m-1,j)+1))
            A[j]=max(A[j],up(min(one,two,comb(m,j)*min(v,NU/m)**j)))
    return A

def reconstruct(n,w):
    coeff=[None]*(n-1)+[Q(0),Q(1)];idx=0
    for _ in range(4):
        for k in range(n-2,1,-1):
            m=n-k;env=caps_by_degree(m);r=Q(w[idx]);idx+=1;E=e_bound(n,k)
            A=(L+E*r)/(2*C0)
            assert r>0
            if r*2*C0<L+E*r:
                total=0
                for q in range(2,m+1):
                    x=env[q]*A**q
                    if coeff[n-q] is not None:x=min(x,coeff[n-q]*Q(comb(m,q),comb(n,q)))
                    total+=x*r**(m-q)
                assert total<r**m,("root witness",n,k,float(total/r**m))
            for q in range(2,m+1):
                c=up(env[q]*A**q*Q(comb(n,q),comb(m,q)))
                i=n-q;coeff[i]=c if coeff[i] is None else min(coeff[i],c)
    assert len(w)==idx
    return coeff

def matrix_variance(n,c,u):
    m=n-1;c=c.copy();c[1]=u
    d=[Q(n-q,n)*c[n-q] for q in range(1,n)]
    # Horner form of M=sum (1-j/n)c_j C^j in the positive companion.
    cur=[[Q(0) for j in range(m)] for i in range(m)]
    for j in range(n-2,-1,-1):
        first=[sum(d[i]*cur[i][h] for i in range(m)) for h in range(m)]
        cur=[first]+cur[:-1]
        if j:
            for i in range(m):cur[i][i]+=Q(n-j,n)*c[j]
    assert all(cur[i][i]>=cur[i+1][i+1] for i in range(m-1))
    off=sum(cur[i][j]*cur[j][i] for i in range(m) for j in range(m) if i!=j)
    diag=sum((cur[i][i]-cur[j][j])**2 for i in range(m) for j in range(i))
    return off+diag/m

def gram_check(n,c,u,R):
    m=n-1;t=root_floor(u/n,m,10**17)
    assert t**m<=u/n and 0<t
    a={j:up(c[j]/t**(n-j)) for j in range(2,n-1)}
    d=[Q(0)]*n;d[m]=Q(1)
    for j in a:d[n-j]=up(Q(j,n)*a[j])
    assert sum(d[q]*R**(m-q) for q in range(1,n))<R**m
    # Positive generating series D'(s)s/(1-D(s)).
    series=[Q(1)]
    for q in range(1,m+1):series.append(up(sum(d[i]*series[q-i] for i in range(1,q+1))))
    powers=[Q(0)]+[up(sum(i*d[i]*series[q-i] for i in range(1,q+1))) for q in range(1,m+1)]
    B=upper_sqrt(sum(q*q for q in d[1:m]))
    need=m*(m-powers[2]);assert need>0
    total=sum(Q(n-j,n)*v*(powers[j-1]+powers[j+1]+B*(1+(R+1)**2*R**(j-1)/4))
              for j,v in a.items())
    return total/need

def circle(n):
    x=Q(1,n);b=(1-x)/2;kappa=L
    v=Q(1);I=Q(1);first=Q(0)
    for j in range(1,161):
        v*=Q(b-j+1,j);I+=v*v
        if j<=10:first+=v*v
    I=up(I+160*v*v)
    H=upper_sqrt((n-1)*I/4);U=2*PHI*(1+(n-1)*I/4)
    # Positive lower split-beta expansion. Replace powers of 2 by outward dyadics.
    lo2=root_floor(Q(2),n);hi2=lo2+Q(1,1<<76)
    w=Q(1);p=Q(1);S=T=Q(0)
    for j in range(55):
        S+=w/(2**j*(2*j+x));T+=p/(2**j*(2*j+1))
        w*=Q(2*j+1,2*j+2);p*=Q(1-x/2+j,j+1)
    J=S/upper_sqrt(hi2)+T/upper_sqrt(Q(2))
    C=2*lo2*J;gap=up(U-C)
    start=Q(1,128*n*n);todo=[(start,Q(1,2)),(Q(1,2),Q(1))];count=0
    while todo:
        l,h=todo.pop()
        gain=0 if h>Q(1,2) else 2*PLO*first*max(Q(0),b*root_floor(1-h*h,2)-H*h)**2
        target=(gap-gain)/kappa
        if target<=0 or target**n<l:count+=1;continue
        assert h-l>Q(1,2**27),("circle",n,float(l),float(h))
        a=l+(h-l)/3;bb=l+2*(h-l)/3
        todo.extend([(l,a),(a,bb),(bb,h)])
    return count

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--low",type=int,default=15);parser.add_argument("--high",type=int,default=27)
    parser.add_argument("--out",default="peeled_independent.json");opt=parser.parse_args();assert __debug__
    l5,h5=atan_inv(5);l239,h239=atan_inv(239)
    assert PLO<16*l5-4*h239<16*h5-4*l239<PHI
    assert Q(433,250)**2<3
    ws=json.loads((G/"peeled_primary.json").read_text())["radii"]
    params=json.loads((G/"collinear_gram_parameters.json").read_text());out={}
    for n in range(opt.low,opt.high+1):
        coeff=reconstruct(n,ws[str(n)]);u,r=map(Q,params[str(n)])
        sm=matrix_variance(n,coeff,u)/2;lg=gram_check(n,coeff,u,r)
        assert sm<Q(94,100) and lg<Q(84,100),("values",n,float(sm),float(lg))
        total=circle(n)
        out[str(n)]=[float(sm),float(lg),total]
        print(n,*out[str(n)],flush=True)
    (G/opt.out).write_text(json.dumps(out,indent=1)+"\n")

