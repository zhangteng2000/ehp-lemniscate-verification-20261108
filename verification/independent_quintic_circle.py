"""Independent, standalone rational audit of the simplified quintic comparison.

It uses coarse scalar bounds, direct double sums for all polynomial pairings,
integer root bisection, and a different bound on the discarded matrix entries.
Only Python's standard library is required.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as F
from math import comb
import json

S=1<<90
def root(q,d):
    q=F(q);assert q>=0 and d>=2
    a=0;b=S
    while F(b**d,S**d)<=q:b*=2
    while b-a>1:
        c=(a+b)//2
        if F(c**d,S**d)<=q:a=c
        else:b=c
    assert F(a**d,S**d)<=q<F((a+1)**d,S**d)
    return F(a,S),F(b,S)

def beta(a,N=70):
    def series(q,z):
        v=F(1);total=F(0)
        for j in range(N):
            total+=v/(z+j);v*=F(q+j,2*(j+1))
        return total,total+F(2,2**N)/(z+N)
    a=F(a);l1,u1=series(F(1,2),a);l2,u2=series(1-a,F(1,2))
    r1,r2=root(2,a.denominator);k=a.numerator
    s1,s2=root(2,2)
    return l1/r2**k+l2/s2,u1/r1**k+u2/s1

def coefficient(q,k):
    v=F(1)
    for j in range(k):v*=F(q+j,j+1)
    return v

def verify():
    pi_lo,pi_hi=F(333,106),F(355,113)
    def atan(q,N):
        v=sum(F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(N))
        return v,v+F(1,(2*N+1)*q**(2*N+1))
    l,u=atan(5,60);a,b=atan(239,16)
    assert 16*l-4*b>pi_lo and 16*u-4*a<pi_hi
    b=F(2,5);q=F(1);moment=F(1);N=200
    for j in range(1,N+1):q*=F(b-j+1,j);moment+=q*q
    mo_lo,mo_hi=F(11831,10000),F(11832,10000)
    assert moment>mo_lo and moment+N*q*q<mo_hi
    assert 2*pi_hi*(1+moment+N*q*q)<F(137169,10000)
    assert root(2,5)[0]*beta(F(1,10))[0]>F(130068,10000)
    gap=F(7101,10000)
    H=F(109,100);H2=mo_hi;assert H*H>H2

    x=F(1,5);R=F(3);tau=F(1,120);s=F(299,100);rc=F(5,3)
    assert s*s*(1+(tau*R)**2)<R*R
    assert beta(F(2,5))[1]<F(369,100)
    assert beta(F(1,10))[1]<F(1133,100)
    M=F(137169,10000)*root(R/2,5)[1]*F(369,100)/pi_lo
    e=tau*R/root(1-(tau*R)**2,2)[0]
    core=M*(e+e*e)+F(6,11)*tau*tau*R*R*root(R,5)[1]
    ratio=(rc/R)**4
    A=H*root(R,5)[1]*(rc/R)**2/(F(9,5)*root(1-ratio,2)[0])
    B=2*H2*root(rc,5)[1]**4/(F(3,5)*root(R,5)[0]**3*(1-ratio))
    motion=(4*A+2*B)*tau/(F(2,5)-tau)
    first=2*F(93,1000)*F(1133,100)*root(tau,5)[1]**4
    assert 1/(8*root(2,2)[0])+F(2,5)*tau<F(93,1000)
    def K(v,l):return 1+1/v+(1+pi_hi*l/2)/min(v-1,1)
    y=1/root(s,5)[0];top=240
    rows=[]
    for j in range(2,top):
        rows.append(sum(K(F(j+k-1,5),F(abs(j-k),5))*y**(j+k)
                        for k in range(2,top) if (j-k)%5==0 and j+k>6))
    m=top+2;g=1/(1-y)
    hsum=y**m*(m*g+y*g*g)
    hsq=y**m*(m*m*g+2*m*y*g*g+y*(1+y)*g**3)
    row=max(rows)+7*hsum+2*hsq
    z=1/s;Q=26
    square=sum((K(F(q)-x,q)*z**q)**2 for q in range(2,Q))
    tail=z**Q*((3+2*Q)/(1-z)+2*z/(1-z)**2)
    A2=root(square,2)[1]+tail
    higher=2*tau*R*root(R,5)[1]*(H2*row+2*H*A2)
    curvature=core+motion+first+higher
    assert curvature<F(9,10)<F(17,16)

    # A pairing is a rational affine function of the single real moment I0.
    def pairing(P,V):
        q=sum(P[i]*V[j]*coefficient(-b,abs(i-j))/coefficient(1+b,abs(i-j))
              for i in range(len(P)) for j in range(len(V)))
        a,c=q*mo_lo-P[0]*V[0],q*mo_hi-P[0]*V[0]
        return min(a,c),max(a,c)
    normals=[];bands=[]
    for m in range(1,4):
        numerator=[coefficient(-b,k)*coefficient(b,m-k) for k in range(m+1)]
        assert sum(numerator)==0
        P=[sum(numerator[:k+1]) for k in range(m)]
        nl,nh=pairing(P,P);assert nl>0;normals.append(nl);band={}
        for J in range(1,m+1):
            for l in range(m-J+1):
                r=m-J-l;V=[F(0)]*m
                for k in range(r+1):
                    a=coefficient(J+l-b,k)*coefficient(b+J,r-k)
                    for j in range(J):V[l+k+j]+=a*comb(J-1,j)*(-1)**j
                a,c=pairing(P,V);cost=max(abs(a),abs(c))
                bound=(F(1,5)+1/(4*root(2,2)[0]) if (J,l)==(1,0) else
                       H*F((J+l)**(J+l),J**J*(l**l if l else 1)))
                band[J+l]=band.get(J+l,F(0))+cost*bound
        bands.append(band)
    ds=[F(480,14401)]+[F(1,k) for k in (24,16,10,6,4,3)]
    margins=[]
    for lo,hi,M in zip(ds,ds[1:],(3,3,3,3,2,1)):
        e=hi/2
        head=sum((1-e*e)**(m+1)*normals[m]-2*sum(v*e**j for j,v in bands[m].items()) for m in range(M))
        assert head>0
        margin=F(17,16)*root(lo/2,5)[0]+2*pi_lo*head-gap
        assert margin>F(3,1000)
        margins.append(F((margin*(1<<40)).__floor__(),1<<40))
    assert F(17,16)*root(F(1,6),5)[0]-gap>F(1,50)
    out={'status':'all standalone rational checks passed',
         'curvature_upper':str(F((curvature*(1<<40)).__ceil__(),1<<40)),
         'curvature_bound':'9/10','minimum_margin_lower':str(min(margins)),
         'moment_interval':[str(mo_lo),str(mo_hi)],'max_rows':3,'intervals':6}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Independent quintic curvature:',float(curvature),'; minimum diameter margin:',float(min(margins)))
    return out

if __name__=='__main__':assert __debug__;verify()
