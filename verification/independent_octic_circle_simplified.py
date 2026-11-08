"""Independent rational degree-eight transfer with two rows and trisection."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from pathlib import Path
import json
from audit_octic import powup,jupper
from audit_peeled_line import PLO,PHI,root_floor,upper_sqrt,atan_inv
from audit_low_degree_transfer import const,head
from audit_octic_majorant_exact import verify as majorant

def verify():
    majorant()
    a,b=atan_inv(5);c,d=atan_inv(239)
    assert PLO<16*a-4*d<16*b-4*c<PHI
    n=8;L=Q(39,50);K=Q(53,100);t=Q(1,300)
    x,b,Il,Ih,gap,H=const(n);R=Q(3);rc=Q(5,3);beta=Q(2,5)
    U=2*PHI*(1+(n-1)*Ih/4)
    e=t*R/root_floor(1-(t*R)**2,2)
    core=U*powup(R/2,x)*2*jupper(1-x)/PLO*(e+e*e)+(1+x)/(2+x)*t*t*powup(R,2+x)
    ratio=(rc/R)**4
    A=H*powup(R,x)*(rc/R)**2/((2-x)*root_floor(1-ratio,2))
    B=2*H**2*powup(rc,4*x)*powup(R,-3*x)/(3*x*(1-ratio))
    motion=(4*A+2*B)*t/(beta-t)
    hn=K/2+Q(7,16)*t
    first=4*hn*jupper(x)*powup(t,1-x)
    s=Q(2999,1000);assert s*s*(1+(t*R)**2)<R*R
    S=1<<44
    def up(v):return Q((v*S).__ceil__(),S)
    y=up(powup(s,-x));cut=208
    def kern(v,k):return 1+1/v+(1+PHI*k/2)/min(1,v-1)
    powers=[Q(1)]
    for _ in range(2*cut+2):powers.append(up(powers[-1]*y))
    h=cut+2;r=1-y
    tail=y**h*((2*h*h+(n+2)*h)/r+(4*h+n+2)*y/r**2+2*y*(1+y)/r**3)
    row=max(sum(kern(Q(i+j-1,n),abs(j-i)//n)*powers[i+j]
                for j in range(2,cut) if (i-j)%n==0 and i+j>n+1)
            for i in range(2,cut))+tail
    z=1/s;q=20
    a2=upper_sqrt(sum((kern(j-x,j)*z**j)**2 for j in range(2,q)))
    a2+=z**q*((3+2*q)/(1-z)+2*z/(1-z)**2)
    higher=2*t*powup(R,1+x)*(H*H*row+2*H*a2)
    curve=core+motion+first+higher
    assert curve<Q(7,10)<L

    data=head(n,2);start=2*t/(1+t*t)
    stack=[(start,Q(3,4)),(Q(3,4),Q(1))];finished=[]
    while stack:
        lo,hi=stack.pop();gain=Q(0)
        if hi<=Q(3,4):
            c=root_floor(1-hi*hi,2);cb=upper_sqrt(1-lo*lo)
            for m,(nl,nh,band) in enumerate(data,1):
                v=max(Q(0),nl*c**m-sum(cost*hi**j*cb**(m-j) for j,cost in band.items()))**2/nh
                if m==1:v=max(v,(Il-1)*max(Q(0),b*c-K*hi)**2)
                gain+=v
        if gap-2*PLO*gain<L*root_floor(lo,n):finished.append((lo,hi));continue
        assert hi-lo>Q(1,1<<24)
        a=(2*lo+hi)/3;z=(lo+2*hi)/3
        stack.extend(((lo,a),(a,z),(z,hi)))
    finished.sort()
    assert finished[0][0]==start and finished[-1][1]==1 and all(a[1]==b[0] for a,b in zip(finished,finished[1:]))
    out={'status':'all independent rational checks passed','tau':'1/300','rows':2,
         'curvature_upper':str(up(curve)),'curvature_bound':'7/10',
         'intervals':[[str(a),str(b)] for a,b in finished],'angular_lower':'Jensen'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Independent degree eight: two rows,',len(finished),'closed intervals; curvature',float(curve))
    return out

if __name__=='__main__':assert __debug__;verify()
