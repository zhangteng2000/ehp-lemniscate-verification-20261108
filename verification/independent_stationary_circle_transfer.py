"""Independent rational trisection and curvature audit in degrees five to eight.

Scalar constants come from separate positive series, roots from integer
bisection, and the angular factor from 101 positive binomial terms.
No primary interval class, diameter partition, or pairing table is used.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')
from pathlib import Path
from fractions import Fraction as Q
import sys,json
from audit_low_degree_transfer import const
from audit_peeled_line import PLO,PHI,root_floor,upper_sqrt,atan_inv
from audit_octic import powup,jupper

def curvature(n,Ih,H,t):
    x=Q(1,n);R=Q(3);r0=Q(5,3);beta=Q(2,5)
    H2=(n-1)*Ih/4;U=2*PHI*(1+H2)
    e=t*R/root_floor(1-(t*R)**2,2)
    core=U*powup(R/2,x)*2*jupper(1-x)/PLO*(e+e*e)+(1+x)/(2+x)*t*t*powup(R,2+x)
    ratio=(r0/R)**4
    A=H*powup(R,x)*(r0/R)**2/((2-x)*root_floor(1-ratio,2))
    B=2*H2*powup(r0,4*x)*powup(R,-3*x)/(3*x*(1-ratio))
    motion=(4*A+2*B)*t/(beta-t)
    k={6:Q(27,100),7:Q(9,25),8:Q(53,100)}
    hn=1/(8*root_floor(2,2))+Q(2,5)*t if n==5 else k[n]/2+Q(n-1,2*n)*t
    first=4*hn*jupper(x)*powup(t,1-x)
    s=Q(2997,1000);assert s*s*(1+(t*R)**2)<R*R
    scale=1<<45
    def up(v):return Q((v*scale).__ceil__(),scale)
    y=up(powup(s,-x));cut=220
    def kern(v,k):return 1+1/v+(1+PHI*k/2)/min(1,v-1)
    powers=[Q(1)]
    for _ in range(2*cut):powers.append(up(powers[-1]*y))
    h=cut+2;r=1-y
    tail=y**h*((2*h*h+(n+2)*h)/r+(4*h+n+2)*y/r**2+2*y*(1+y)/r**3)
    row=max(sum(kern(Q(i+j-1,n),abs(j-i)//n)*powers[i+j]
                for j in range(2,cut) if (i-j)%n==0 and i+j>n+1)
            for i in range(2,cut))+tail
    z=1/s;q=24
    a2=upper_sqrt(sum((kern(j-x,j)*z**j)**2 for j in range(2,q)))
    a2+=z**q*((3+2*q)/(1-z)+2*z/(1-z)**2)
    higher=2*t*powup(R,1+x)*(H2*row+2*H*a2)
    return core+motion+first+higher

def verify(n):
    x,b,Il,Ih,gap,H=const(n)
    H2=(n-1)*Ih/4
    K={5:Q(9,50),6:Q(27,100),7:Q(9,25),8:Q(53,100)}[n]
    if n==5:assert K*K*32>1
    L={5:Q(17,16),6:Q(73,100),7:Q(367,500),8:Q(39,50)}[n]
    tau=Q(1,{5:120,6:200,7:300,8:300}[n]);start=2*tau/(1+tau*tau)
    angle=Q(1);coef=Q(1)
    for j in range(1,101):coef*=Q(x/2+j-1,j);angle+=coef*coef
    loss={5:Q(9,10),6:Q(681,1000),7:Q(131,250),8:Q(7,10)}[n]
    curve=curvature(n,Ih,H,tau);assert curve<loss<L
    first=Il-1+b*b/H2
    t=b*b*(1-b)*(b+3)/4;s=b*(b+1)/2
    edge=Q(173,200);todo=[(start,edge),(edge,Q(1))];done=[]
    while todo:
        lo,hi=todo.pop();cost=Q(0)
        if hi<=edge:
            c=root_floor(1-hi*hi,2);cb=upper_sqrt(1-lo*lo)
            g=1/root_floor(2*(1+c),2)
            hn=max(0,b*c-K*g*hi)
            d=max(0,t*(c*c-2*K*g*hi*cb/b)-hi*hi*(b*((n-1)*x/2+8*K/3)+s*H))
            cost=2*PLO*(first*hn*hn+d*d)
        if cost+L*angle*root_floor(lo,n)>gap:done.append((lo,hi));continue
        assert hi-lo>Q(1,2**27),(n,str(lo),str(hi))
        a=(2*lo+hi)/3;z=(lo+2*hi)/3;todo.extend(((lo,a),(a,z),(z,hi)))
    done.sort();assert done[0][0]==start and done[-1][1]==1 and all(a[1]==b[0] for a,b in zip(done,done[1:]))
    print('Independent exact stationary transfer',n,len(done),'closed intervals',flush=True)
    return {'n':n,'nonlinear_bound_used':str(K),'angular_sum_terms':101,'curvature_upper':str(curve),'curvature_bound':str(loss),
            'intervals':[[str(a),str(b)] for a,b in done]}

if __name__=='__main__':
    a,b=atan_inv(5);c,d=atan_inv(239)
    assert PLO<16*a-4*d<16*b-4*c<PHI
    out={str(n):verify(n) for n in range(5,9)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
