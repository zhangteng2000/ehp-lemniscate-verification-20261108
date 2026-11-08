"""Rational orthogonal-projection bounds for the all-degree collinear reduction."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as F
from pathlib import Path
import json
import argparse
from audit_quintic_rational import I,PI,J,taylor,vin,iv

def constants(n):
    x=F(1,n);b=(1-x)/2
    C=I(2)**(1+x)*J(x);M=I(2)**(2-x)*J(2-x)/PI
    U=2*PI*(1+(1-x)/(4*x)*M)
    # angular lower multiplier for a farthest pair, minimized at antipodes
    Jmin=I(2)**(1-x)*J(1-x)/PI
    return x,b,C,M,U,Jmin

def curvature(n,R,tau,prop):
    x,b,C,m,U,jmin=prop
    R=I(R);tau=I(tau);z=tau*R;s=R/(1+z*z)**F(1,2)
    eps=z/(1-z*z)**F(1,2);M=U*(R/2)**x*2*J(1-x)/PI
    core=M*(eps+eps*eps)+(1+x)/(2+x)*tau**2*R**(2+x)
    straight=2*n*R**(-2+x)/((2*n-1)*(1-R**-4))**F(1,2)
    v=s**-2;Y=((2*n-1)/(1-v)+n*v/(1-v)**2)**F(1,2)
    W=max((s**(-2*(j-1)*x)*(j+(n-j)*v)/(1-v)**2).h for j in range(1,n+1))
    curved=I(n-1)**F(1,2)*2*J(x)*tau**(1-x)+2*Y*s**(-2+x)+W*s**(-3*x)/min(3*x,1)
    return (core+straight+curved).h

def rows(n,M,constants):
    x,b,C,moment,U,jmin=constants
    H=((1-x)/(4*x)*moment)**F(1,2)
    K={};out=[]
    for m in range(1,M+1):
        P=[];v=F(0)
        for k in range(m):
            v+=taylor(-b,k)*taylor(b,m-k);P.append(v)
        assert v+taylor(-b,m)==0
        p=[sum(P[i]*taylor(-b,abs(i-j))/taylor(1+b,abs(i-j)) for i in range(m)) for j in range(m)]
        def pair(v):return moment*sum(a*b for a,b in zip(p,v))-P[0]*v[0]
        N=pair(P);assert N.l>0;band={}
        for j in range(1,m+1):
            for l in range(m-j+1):
                if (j,l) not in K:
                    val=H.h*F((j+l)**(j+l),j**j*(l**l if l else 1))
                    if (j,l)==(1,0):
                        if n==5:val=min(val,(F(1,5)+1/(4*I(2)**F(1,2))).h)
                        if n==6:val=min(val,F(21,40))
                    K[j,l]=val
                cost=pair(vin(m,j,l,b)).absupper()*K[j,l]
                band[j+l]=band.get(j+l,0)+cost
        # Project each unknown polynomial onto the model vector S_m and
        # retain the square of the surviving projection. This bounds
        # every row from below by a nonnegative quantity.
        out.append((N,band))
    return out

def lower(lo,hi,line,n,data,const):
    x,b,C,moment,U,jmin=const
    tail=I(0);used=0
    modlo=(1-(hi/2)**2)**1 if hi<=F(3,2) else 0
    # The chord disk is only needed when hi<=3/2<sqrt(3).
    if hi<=F(3,2):
        modhi=I(1-(lo/2)**2)**F(1,2);modlo=I(modlo)**F(1,2)
        for m,(N,coef) in enumerate(data,1):
            correction=sum(val*(hi/2)**j*modhi.h**(m-j) for j,val in coef.items())
            retained=max(F(0),modlo.l**m*N.l-correction)
            if retained>0:used+=1
            tail+=retained**2/N.h
    return (line*I(lo/2)**x*jmin+2*PI*tail-(U-C)).l,used

def run(n,M,tau,R,line,diagonal=False):
    const=constants(n);data=rows(n,M,const)
    if diagonal:
        from diagonal_curvature_bound import certified
        k=certified(n,R,tau,F(2,5),const)
    else:k=curvature(n,R,tau,const)
    assert k<line,("curvature",float(k),float(line))
    d0=4*tau/(1+tau*tau)
    todo=[(d0,F(2),0)];proved=[];least=None
    while todo:
        lo,hi,depth=todo.pop();v,used=lower(lo,hi,line,n,data,const)
        if v>F(1,100000):
            proved.append([str(lo),str(hi),used]);least=v if least is None else min(v,least)
        else:
            if depth>=19:
                raise ValueError(("open interval",float(lo),float(hi),float(v),used))
            mid=(lo+hi)/2;todo.extend([(mid,hi,depth+1),(lo,mid,depth+1)])
    out={"n":n,"M":M,"R":R,"tau":str(tau),"diagonal_curvature":diagonal,"diameter_two_line":str(line),
         "curvature_upper_display":float(k),"minimum_margin_gt":"1/100000","minimum_margin_display":float(least),
         "closed_intervals":proved}
    fil=Path(__file__).with_name(f"coherent_projection_{n}_{line.numerator}_{line.denominator}.json")
    fil.write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
    print(n,"line",line,"small curvature",float(k),"diameter intervals",len(proved),"minimum margin",float(least))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--n",type=int,default=6);p.add_argument("--M",type=int,default=19)
    p.add_argument("--tau",type=F,default=F(1,3000));p.add_argument("--R",type=int,default=25)
    p.add_argument("--diagonal",action="store_true")
    p.add_argument("--line",type=F,default=F(4,5));a=p.parse_args();run(a.n,a.M,a.tau,a.R,a.line,a.diagonal)

