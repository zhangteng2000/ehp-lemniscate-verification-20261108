"""Exact two-coefficient circle transfer in degrees five through eight.

The analytic first-variation orthogonality and affine inverse identity
are proved in the paper. This program certifies their scalar consequences
on closed intervals and checks the common short-arc curvature estimate.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')
from pathlib import Path
from fractions import Fraction as Q
import json,sys
from coherent_projection_transfer import constants
from audit_quintic_rational import I,PI
from diagonal_curvature_bound import certified

def run(n):
    cc=constants(n);x,b,C,I0,U,J=cc
    L={5:Q(17,16),6:Q(73,100),7:Q(367,500),8:Q(39,50)}[n]
    K=(1/(4*I(2)**Q(1,2))).h if n==5 else {6:Q(27,100),7:Q(9,25),8:Q(53,100)}[n]
    tau=Q(1,{5:120,6:200,7:300,8:300}[n]);start=2*tau/(1+tau*tau)
    loss={5:Q(9,10),6:Q(681,1000),7:Q(131,250),8:Q(7,10)}[n]
    hn=(1/(8*I(2)**Q(1,2))).h+Q(2,5)*tau if n==5 else K/2+Q(n-1,2*n)*tau
    curve=certified(n,3,tau,Q(2,5),cc,hn_bound=hn)
    assert curve<loss<L
    H2=(U/(2*PI)-1).h;H=(I(H2)**Q(1,2)).h
    factor=I0.l-1+b*b/H2
    T=-b*b*(b-1)*(b+3)/4;S=b*(b+1)/2
    edge=Q(173,200);todo=[(start,edge),(edge,Q(1))];done=[];least=None
    while todo:
        lo,hi=todo.pop();gain=Q(0)
        if hi<=edge:
            c=(I(1-hi*hi)**Q(1,2)).l;cb=(I(1-lo*lo)**Q(1,2)).h
            g=(1/(I(2*(1+c))**Q(1,2))).h
            h=max(0,b*c-K*g*hi)
            z=max(0,T*(c*c-2*K*g*hi*cb/b)-hi*hi*(b*(Q(n-1,2*n)+8*K/3)+S*H))
            gain=2*PI.l*(factor*h*h+z*z)
        target=((U-C).h-gain)/(L*J.l)
        if target<=0 or target**n<lo:
            done.append((lo,hi));continue
        assert hi-lo>Q(1,2**26),(n,str(lo),str(hi))
        mid=(lo+hi)/2;todo.extend(((lo,mid),(mid,hi)))
    done.sort();assert done[0][0]==start and done[-1][1]==1 and all(a[1]==b[0] for a,b in zip(done,done[1:]))
    out={'n':n,'tau':str(tau),'line_constant':str(L),'nonlinear_cap_upper':str(K),
         'curvature_upper':str(Q((curve*(1<<60)).__ceil__(),1<<60)),'curvature_bound':str(loss),
         'intervals':[[str(a),str(b)] for a,b in done]}
    print('Exact stationary transfer',n,len(done),'closed intervals',flush=True)
    return out

if __name__=='__main__':
    out={str(n):run(n) for n in range(5,9)}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
