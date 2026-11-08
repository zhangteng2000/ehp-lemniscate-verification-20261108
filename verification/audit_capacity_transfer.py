"""Independent arithmetic for the circular-arc capacity transfer, including the cubic."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
from functools import lru_cache
from math import comb
from audit_quintic_rational import I,PI,taylor,vin
from coherent_projection_transfer import constants
from audit_peeled_line import PLO,PHI,root_floor,upper_sqrt
from audit_octic import powup,jupper
import json,argparse
G=Path(__file__).resolve().parent;SS=1<<45
def up(y):return Q((SS*y).__ceil__(),SS)
def dn(y):return -up(-y)
def area_ball(n,goal):
    if n==3:return 0
    cells=[((Q(0),Q(0)),(Q(1),Q(1)))];num=0
    while cells:
        a,b=cells.pop();U,V=b;u,v=a
        if u+v>1:continue
        if n==6:
            z=upper_sqrt(3*U*(1-u-v))/4+V/4+5*upper_sqrt(U**3)/48
        else:
            z=upper_sqrt((4*U+6*V)*(1-u-v))/4+5*upper_sqrt(2*V)*U/16
        if z<goal:num+=1;continue
        k=0 if U-u>V-v else 1;w=(a[k]+b[k])/2
        assert b[k]-a[k]>Q(1,16384);l=list(a);r=list(b);l[k]=w;r[k]=w
        cells.extend([(tuple(l),b),(a,tuple(r))])
    return num
def curvature(n,L,t,c,signed=False,hn_bound=None):
    x=Q(1,n);H=upper_sqrt((n-1)*c[3].h/4);U=c[4].h;R=Q(3);r0=Q(5,3);b=Q(2,5)
    ep=t*R/root_floor(1-(t*R)**2,2)
    core=U*powup(R/2,x)*2*jupper(1-x)/PLO*(ep+ep*ep)+(1+x)/(2+x)*t*t*powup(R,2+x)
    rr=(r0/R)**4
    move=(4*H*powup(R,x)*(r0/R)**2/((2-x)*root_floor(1-rr,2))
          +4*H*H*powup(r0,4*x)*powup(R,-3*x)/(3*x*(1-rr)))*t/(b-t)
    hn=min(H,Q(21,40)+t/6) if n==6 else H
    if signed:hn=min(hn,{6:Q(27,200),7:Q(18,100)}[n]+Q(n-1,2*n)*t)
    if hn_bound is not None:hn=min(hn,Q(hn_bound))
    first=4*hn*jupper(x)*powup(t,1-x)
    s=Q(2997,1000);assert s*s*(1+(t*R)**2)<R*R
    y=up(powup(s,-x));N=220;v=1/s
    def k(a,z):return 1+1/a+(1+PHI*z/2)/min(1,a-1)
    yy=[Q(1)]
    for _ in range(2*N):yy.append(up(y*yy[-1]))
    row=max(sum(k(Q(j+i-1,n),abs(i-j)//n)*yy[i+j] for i in range(2,N)
                if (i-j)%n==0 and i+j>n+1) for j in range(2,N))
    h=N+2;p=1-y;row+=y**h*((2*h*h+(n+2)*h)/p+(4*h+n+2)*y/p**2+2*y*(1+y)/p**3)
    J=24
    a2=upper_sqrt(sum((k(i-x,i)*v**i)**2 for i in range(2,J)))+v**J*((3+2*J)/(1-v)+2*v/(1-v)**2)
    high=2*t*powup(R,1+x)*(H**2*row+2*H*a2);assert core+move+first+high<L
    return float(core+move+first+high)
def form(n,M,c):
    x,b,_,I0,_,_=c;H=upper_sqrt((n-1)*I0.h/4);out=[]
    for m in range(1,M+1):
        P=[];v=Q(0)
        for j in range(m):v+=taylor(-b,j)*taylor(b,m-j);P.append(v)
        T=[sum(P[i]*taylor(-b,abs(i-j))/taylor(1+b,abs(i-j)) for i in range(m)) for j in range(m)]
        def ip(W):
            z=I0*sum(p*h for p,h in zip(T,W))-P[0]*W[0]
            return I(dn(z.l),up(z.h))
        N=ip(P);assert N.l>0;one=[ip(vin(m,1,i,b)) for i in range(m)];other=[Q(0)]*(m+1)
        for j in range(2,m+1):
            for i in range(0,m+1-j):
                a=ip(vin(m,j,i,b));K=H*Q((j+i)**(j+i),j**j*(i**i if i else 1))
                other[j+i]+=up(a.absupper()*K)
        out.append((N,one,other))
    return out
@lru_cache(None)
def coefficients(x):
    lam=x/2
    return [[taylor(lam,i)*taylor(lam,k-i) for i in range(k+1)] for k in range(29)]
def chord(x,c0,c1):
    mid=(c0+c1)/2;rad=(c1-c0)/2;cos=[Q(1),mid]
    for k in range(2,29):cos.append(2*mid*cos[-1]-cos[-2])
    s=Q(0)
    for k,a in enumerate(coefficients(x)):
        val=sum(w*cos[abs(k-2*i)] for i,w in enumerate(a))
        err=rad*sum(w*(k-2*i)**2 for i,w in enumerate(a))
        s+=dn(max(0,abs(val)-err)**2)
    return s
def circle(n,L,M,t,goal,c):
    x,b,C,I0,U,J=c;data=form(n,M,c);edge=Q(173,200);start=2*t/(1+t*t)
    stack=[(start,edge),(edge,Q(1))];cnt=0;gap=(U-C).h;rho=[Q(j,40) for j in range(1,40)]
    while stack:
        lo,hi=stack.pop();z=dn(root_floor(1-hi*hi,2));Z=min(1,up(upper_sqrt(1-lo*lo)))
        gain=Q(0);jmin=max(J.l,chord(x,z,Z))
        if hi<=edge:
            g=up(1/root_floor(2*(1+z),2));p=hi/(1+z);br=[]
            for l in range(M):
                if l==0:d=goal*g
                else:d=goal*min(min(1/(1-v*v),g/(1-v)**2)/v**l for v in rho)
                R=Q(n-3,2*n)
                if l%2:
                    endpoint=Q((-1)**((l+1)//2),n);a=endpoint-R;e=endpoint+R
                elif l==0:a=0;e=R*p
                else:e=R*min(1,(l+1)*p);a=-e
                br.append(I(dn(a-d),up(e+d)))
            for m,(N,one,other) in enumerate(data,1):
                lower=dn(N.l*z**m)
                for i,a in enumerate(one):
                    t0=I(lo**(i+1)*z**(m-1-i),hi**(i+1)*Z**(m-1-i))*br[i]*a
                    lower+=dn(t0.l)
                lower-=sum(up(b*hi**j*Z**(m-j)) for j,b in enumerate(other))
                gain+=dn(max(0,lower)**2/N.h)
        if gap<2*PI.l*gain+L*jmin*root_floor(lo,n):cnt+=1;continue
        assert hi-lo>Q(1,1<<26),(n,float(lo))
        b1=(2*lo+hi)/3;b2=(lo+2*hi)/3;stack.extend([(lo,b1),(b1,b2),(b2,hi)])
    return cnt
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,nargs='+',default=[6,7]);p.add_argument('--simplified',action='store_true');a=p.parse_args()
    assert __debug__;pars={3:(Q(6,5),15,Q(1,110),Q(0)),6:(Q(73,100),20,Q(1,230),Q(27,100)),7:(Q(367,500),16,Q(1,392),Q(36,100))}
    if a.simplified:
        pars[6]=(Q(73,100),12,Q(1,200),Q(27,100))
        pars[7]=(Q(367,500),6,Q(1,300),Q(36,100))
    out={}
    for n in a.n:
        L,M,t,goal=pars[n];con=constants(n)
        if a.simplified:
            from audit_capacity_majorants_exact import verify as exact_majorants
            exact_majorants();v=0
        else:v=area_ball(n,goal)
        w=curvature(n,L,t,con,signed=a.simplified);k=circle(n,L,M,t,goal,con)
        out[n]=[v,w,k];print(n,'unit-ball leaves',v,'curvature',w,'diameter intervals',k,flush=True)
    (G/('capacity_transfer_simplified_independent.json' if a.simplified else 'capacity_transfer_independent.json')).write_text(json.dumps(out,indent=1)+'\n')

