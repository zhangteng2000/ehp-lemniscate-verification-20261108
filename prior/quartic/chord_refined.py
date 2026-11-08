"""Implicit analytic enclosures at the two kinds of collision and a signed exterior bound."""
from fractions import Fraction as F
import math
import slice_data as sd
from chord_geometry import (J, arb, acb, I, OMEGA, ALPHA, F16, real, interval, rational, fixed_down,
                            fixed_up, algebra, Outside, Split)

TRUNC=16
SQRT3=arb(3).sqrt()


def entry(ss,tt,hs,ht):
    sr,tr=real(ss),real(tt)
    t=interval(tt,ht)
    h0=(1+I*tr)/(1-I*tr)
    # The complex unit h stays on the circle along the real parameter interval.
    h=h0+acb(arb(0,2*real(ht)),arb(0,2*real(ht)))
    wh0=h0**3
    sf=(6*tr.atan()).sin();cf=(6*tr.atan()).cos()
    angle=6*t.atan()
    sn=angle.sin();cs=angle.cos()
    Y=J(sn,0,6*cs/(1+t*t));X=J(cs,0,-6*sn/(1+t*t))
    r=interval(ss,hs);rho=J(r*r*r,3*r*r,0)
    v=rho*(X+I*Y)
    u=J(r*h,h,2*I*(r*h)/(1+t*t))
    return u,rho,X,Y,v,sr*h0,sr**3*wh0


def certify_disc(center,lin,higher,residual,variable,perturb):
    # Scalar Rouche bound in the unknown for all input parameters in one disk.
    slope=abs(lin)
    if not slope>0:raise Split("implicit derivative")
    need=residual+variable*perturb[0]
    radius=(2*need/slope).upper()+arb(2)**-96
    for _ in range(4):
        rhs=residual+sum(abs(c)*radius**degree for degree,c in enumerate(higher,2))
        rhs+=variable*sum(v*radius**degree for degree,v in enumerate(perturb))
        if slope*radius>rhs:
            ball=arb(0,radius)
            return center+acb(ball,ball)
        radius*=2
    raise Split("implicit root disk")


def shape_root(sh,ss,tt,hs,ht,sv,u,v,u0,v0,c0):
    ds,dt=real(hs),real(ht)
    if sh==0:
        dv=(real(ss)+ds)**3-real(ss)**3+6*real(ss)**3*dt
        p0=c0**4-6*c0**2+(16*v0-8)*c0-3
        lin=4*c0**3-12*c0+16*v0-8
        c=certify_disc(c0,lin,[6*c0**2-6,4*c0,1],abs(p0),dv,[16*abs(c0),16])
        denom=4*c**3-12*c+16*v.v-8
        if not abs(denom)>0:raise Split("implicit simple root")
        cv=-16*c/denom
        return J(c,cv*v.s,cv*v.t)
    if ss:
        q0=(c0+1)/u0
    else:q0=-arb(4).root(3)*[0,1,OMEGA,OMEGA.conjugate()][sh]
    du=ds+2*real(ss)*dt
    p0=4*q0**3-u0*q0**4+16-16*u0*q0
    lin=12*q0**2-4*u0*q0**3-16*u0
    q=certify_disc(q0,lin,[12*q0-6*u0*q0**2,4-4*u0*q0,-u0],abs(p0),du,
                   [abs(q0**4+16*q0),abs(4*q0**3+16),6*abs(q0)**2,4*abs(q0),1])
    denom=12*q**2-4*u.v*q**3-16*u.v
    if not abs(denom)>0:raise Split("implicit blown-up root")
    qu=(q**4+16*q)/denom
    cu=q+u.v*qu
    return J(-1+u.v*q,cu*u.s,cu*u.t)


def vectors(c,X,Y,r,full):
    m=c*F16*(Y+(r-X)*I)*(1/ALPHA**4)
    if full:m=m/Y; lead=J(1)
    else:lead=Y
    mc=m.conj()*c;norm=c.norm2()
    return [lead,2*lead*c.re(),-2*lead*c.im(),lead*norm,-m.re(),-m.im(),-mc.re(),mc.im()]


def mag_pair(j,j0,ds,dt):
    shift=abs(j.s)*ds+abs(j.t)*dt
    a=abs(j.v);b=abs(j0.v)
    return (max(F(0),rational(a.lower()),rational((b-shift).lower())),
            min(rational(a.upper()),rational((b+shift).upper())))


def simple_clip(c,c0,ds,dt):
    norm=c.norm2();n0=c0.norm2()
    lq=[J(1),2*c.re(),-2*c.im(),norm, -norm*(TRUNC**2*F16**2/ALPHA**8)]
    lp=[J(1),2*c0.re(),-2*c0.im(),n0, -n0*(TRUNC**2*F16**2/ALPHA**8)]
    cc=[];cr=[]
    for p,f in zip(lp,lq):
        err=abs(f.s.real)*ds+abs(f.t.real)*dt
        fac=1<<12
        l=max(fixed_down(f.v.real*fac),fixed_down((p.v.real-err)*fac))
        h=min(fixed_up(f.v.real*fac),fixed_up((p.v.real+err)*fac))
        mid=(l+h)//2;cc.append(mid);cr.append(max(mid-l,h-mid))
    return cc,cr


def geometry(sh,ss,tt,hs,ht):
    ds,dt=real(hs),real(ht)
    u,r,X,Y,v,u0,v0=entry(ss,tt,hs,ht)
    if v.v.real>arb(1)/2 or X.v.real<0 or r.v.real>2*X.v.real:raise Outside
    sr,tr=real(ss),real(tt)
    X0=J((6*tr.atan()).cos());Y0=J((6*tr.atan()).sin());r0=J(sr**3);sv0=J(sr)
    vv0=r0*(X0+I*Y0)
    c0=algebra(sh,sv0,J(tr))[2]
    c=shape_root(sh,ss,tt,hs,ht,J(interval(ss,hs),1,0),u,v,u0,v0,c0.v)
    cl,cu=mag_pair(c,c0,ds,dt);cl=max(cl,F(1,6))
    dl,du=mag_pair(1-v,1-vv0,ds,dt);dl=max(dl,F(1,2));du=min(du,F(1))
    if dl>du:raise Outside
    kl=max(F(0),2*rational(Y.v.real.lower())/du)
    kh=2*min(F(1),rational(Y.v.real.upper()))/dl
    if kh<kl:raise Outside
    full=kl>=F(1,12)
    cv=vectors(c,X,Y,r,full);center=vectors(c0,X0,Y0,r0,full)
    model_lower=F(11069,1000)
    if full:target=model_lower*sd.fourthlo(1/kh)
    else:
        q=min(F(1),TRUNC*kl/2);q=F((q*(1<<30)).__floor__(),1<<30)
        val=F(0)
        for a in reversed(sd.TWEIGHTS):val=val*q*q+a
        target=16*F((val*(1<<38)).__floor__(),1<<38)
        if kh*TRUNC>=2:target=min(target,model_lower*sd.fourthlo(1/kh))
        # |a_3| in the monic exterior inverse. The normalizing fourth-root phase cancels.
        c2=c.v*c.v
        U=-2-2*c2/3
        C=1+8*c.v/3+10*c2/9-c2*c2/27
        if abs(c.v)>0:
            a3max=min(F(3,5),rational(abs((-C/4+U*U/32)/(F16*c.v)).upper()))
        else:a3max=F(3,5)
        powk=sd.sqrtup(kh)*sd.sqrtup(sd.sqrtup(kh))
        lead=a3max*min(F(1),18*powk)
        tail=(F(13,25)+F(1,6))/3 + lead
        target-=tail
    target*=sd.fourthlo(F(16,3)*cl)*F(2,5)
    norm=max(float(abs(row.v).upper()) for row in center)
    if not math.isfinite(norm) or norm<=0:raise Split("scale")
    _,e=math.frexp(norm);scale=arb(2)**(sd.SHAPE_BITS-e)
    gens=[[],[],[]];rad=[]
    for p,f in zip(center,cv):
        terms=[p.v.real*scale,f.s.real*ds*scale,f.t.real*dt*scale]
        err=arb(0)
        for g,t in zip(gens,terms):
            a=fixed_down(t.mid());g.append(a);err+=abs(t-a)
        rad.append(fixed_up(err))
    cc,cr=simple_clip(c,c0,ds,dt)
    if sum(sum(map(abs,row)) for row in gens)+sum(rad)>1<<32:raise Split("coefficient width")
    if sum(map(abs,cc))+sum(cr)>1<<32:raise Split("clip width")
    if full:
        O=(r-X)/Y
        scal=1+abs(O.v)
        scs=ds*(abs(c.s)/(1+cl.numerator/cl.denominator)+abs(O.s)/scal)
        sct=dt*(abs(c.t)/(1+cl.numerator/cl.denominator)+abs(O.t)/scal)
    else:
        # A change in inverse target radius changes core angles by O(TRUNC delta k).
        scs=ds*(abs(c.s)+TRUNC*abs(v.s)*abs(Y.v))
        sct=dt*(abs(c.t)+TRUNC*(abs(Y.t)+abs(v.t)*abs(Y.v)))
    scores=float(scs.upper()),float(sct.upper())
    return full,F((target*(1<<38)).__floor__(),1<<38),gens,rad,cc,cr,scores,kl,kh


