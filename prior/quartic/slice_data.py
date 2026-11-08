"""Exact source controls, quartic shape charts, and comparison constants."""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, lcm, isqrt
from pathlib import Path
import argparse
import json
import math
import struct
import numpy as np

ROOT_BITS = 36
SHAPE_BITS = 25
HALF = 4
T = 81
MODEL_LOWER = F(551, 50)
TWEIGHTS = tuple(F(comb(2*m,m), 4**m*(8*m+1)) for m in range(25))


def directions(count=8):
    # Each listed (x,y,d) is a rational unit vector (x+i*y)/d.
    if count == 8:
        first=[(1,0,1),(12,5,13),(119,120,169),(5,12,13)]
    elif count == 12:
        first=[(1,0,1),(60,11,61),(15,8,17),(119,120,169),(8,15,17),(11,60,61)]
    elif count == 4:
        first=[(1,0,1),(119,120,169)]
    else:
        raise ValueError("unsupported polygon")
    grid=first+[(-y,x,h) for x,y,h in first]
    return [(264*x-23*y,23*x+264*y,265*h) for x,y,h in grid]


def weights(dirs):
    gaps=[]
    x0,y0,q0=dirs[0]
    for (x,y,p),(z,w,q) in zip(dirs,dirs[1:]+[(-x0,-y0,q0)]):
        assert x*x+y*y==p*p and x*w-y*z>0
        gaps.append(F(x*w-y*z,p*q+x*z+y*w))
    return [(gaps[i-1]+v)/2 for i,v in enumerate(gaps)]


def add(a,b):
    c=dict(a)
    for i,v in b.items():
        u=c.get(i,(0,0));c[i]=(u[0]+v[0],u[1]+v[1])
    return {k:v for k,v in c.items() if v!=(0,0)}


def mul(a,b):
    c={}
    for (i,j),(u,v) in a.items():
        for (k,h),(x,y) in b.items():
            ij=i+k,j+h
            old=c.get(ij,(0,0))
            c[ij]=(old[0]+u*x-v*y,old[1]+u*y+v*x)
    return {k:v for k,v in c.items() if v!=(0,0)}


def sc(a,k):
    return {i:(x*k,y*k) for i,(x,y) in a.items() if k}


def cj(a):
    return {i:(x,-y) for i,(x,y) in a.items()}


def re(a):
    return {i:(x,0) for i,(x,y) in a.items() if x}


def im(a):
    return {i:(y,0) for i,(x,y) in a.items() if y}


def ip(a):
    return {i:(-y,x) for i,(x,y) in a.items()}


def cn(t):
    return {(0,0):(t,0)} if t else {}


def powp(a,n):
    b=cn(1)
    for _ in range(n):b=mul(a,b)
    return b


def matrix(poly,axis_degree,part=0):
    out=np.zeros((axis_degree+1,axis_degree+1),dtype=object)
    for key,z in poly.items():out[key]=z[part]
    return out


def source(d, family="b"):
    x,y,h=d
    z={(1,0):(F(x,h),F(y,h)),(0,1):(F(-y,h),F(x,h))}
    if family=="chord":
        # The physical source coordinate is (5/2)Z, and both factors
        # have been divided by (5/2)^4.
        a=F(2,5)
        q0=powp(add(powp(z,2),cn(-a*a)),2)
        q1=sc(mul(add(z,cn(-2*a)),powp(add(z,cn(a)),2)),-4*a/3)
    else:
        q0=powp(add(powp(z,2),cn(-1)),2)
        q1=powp(add(z,cn(-1)),2)
    cross=mul(q1,cj(q0))
    return [matrix(p,8,k) for p,k in [(mul(q0,cj(q0)),0),(cross,0),(cross,1),
                                      (mul(q1,cj(q1)),0),(q0,0),(q0,1),(q1,0),(q1,1), (cn(1),0)]]


def build_source(path,nx,ndir,family="b"):
    dirs=directions(ndir)
    # X=4(2j+1-nx)/nx. The Y blossom has 8-i copies of -4, i copies of 4.
    xx=np.array([[((2*j+1-nx)*HALF)**p*nx**(8-p) for p in range(9)]
                 for j in range(nx)],dtype=object)
    m=lcm(*(comb(8,p) for p in range(9)))
    ys=np.array([[sum(comb(i,k)*comb(8-i,p-k)*(-HALF)**(p-k)*HALF**k
                      for k in range(max(0,p-8+i),min(i,p)+1))*(m//comb(8,p))
                  for i in range(9)] for p in range(9)],dtype=object)
    with Path(path).open("wb") as out:
        out.write(struct.pack("<4i",nx,ndir,ROOT_BITS,HALF))
        for d in dirs:
            controls=[]
            for mat in source(d, family):
                den=lcm(*(F(v).denominator for v in mat.flat))
                integers=np.array([int(x*den) for x in mat.flat],dtype=object).reshape(9,9)
                raw=((xx @ integers @ ys)*(1<<ROOT_BITS))//(nx**8*den*m)
                assert max(abs(q) for q in raw.flat)<2**60
                controls.append(np.array(raw,dtype="<i8"))
            np.stack(controls,axis=1).tofile(out)


@lru_cache(None)
def shape_polys(chart):
    # Four charts remove the homogeneous zeros at b=0 and -1/2.
    if chart=="raw":
        b={(1,0):(1,0),(0,1):(0,1)};order=0
    else:
        left=chart[0]=="L"
        sign=1 if left else -1
        if chart[1]=="H": b={(1,0):(sign,0),(1,1):(0,1)}
        else:b={(1,1):(sign,0),(1,0):(0,1)}
        if left:b=add(b,cn(F(-1,2)))
        order=3 if left else 1
    be=mul(b,add(cn(1),b))
    d=mul(powp(add(cn(1),sc(b,-1)),3),add(cn(1),b))
    e=sc(mul(b,powp(add(b,cn(2)),3)),-1)
    f=sc(powp(add(sc(b,2),cn(1)),3),-1)
    det=im(mul(e,cj(d)));on=ip(mul(mul(d,e),cj(f)))
    prod=mul(cj(on),be);be2=re(mul(be,cj(be)))
    pol=[sc(det,2),sc(mul(det,re(be)),-8),sc(mul(det,im(be)),8),
         sc(mul(det,be2),8),sc(re(on),-2),sc(im(on),-2),sc(re(prod),4),sc(im(prod),-4)]
    clean=[]
    for p in pol:
        assert all(i>=order for i,j in p)
        clean.append({(i-order,j):v for (i,j),v in p.items()})
    clean.extend([re(be),im(be),be2])
    g=max(max(i,j) for p in clean for i,j in p)
    den=lcm(*(F(a).denominator for p in clean for a,y in p.values()))
    return g,den,[matrix({key:(int(a*den),0) for key,(a,y) in p.items()},g) for p in clean]


def translate(c,h,g,den):
    x=int(c*den);width=int(h*den)
    return np.array([[comb(p,i)*x**(p-i)*width**i*den**(g-p) if p>=i else 0
                      for p in range(g+1)] for i in range(g+1)],dtype=object)


def sqrtlo(x,bits=44):
    a=1<<bits
    return F(isqrt((x.numerator*a*a)//x.denominator),a)


def sqrtup(x,bits=44):
    a=sqrtlo(x,bits)
    return a if a*a==x else a+F(1,1<<bits)


def fourthlo(v):
    return sqrtlo(sqrtlo(v))


def abs_range(center,rad):
    return max(0,abs(center)-rad),abs(center)+rad


def compact_shape(chart,x,y,hx,hy,order=1):
    g,den,pol=shape_polys(chart)
    H=lcm(x.denominator,y.denominator,hx.denominator,hy.denominator)
    xx=translate(x,hx,g,H); yy=translate(y,hy,g,H).T
    coeff=np.array([xx@mat@yy for mat in pol],dtype=object)
    total=abs(coeff).sum(axis=(1,2))
    al,ah=abs_range(coeff[0,0,0],total[0]-abs(coeff[0,0,0]))
    pl,ph=abs_range(coeff[4,0,0],total[4]-abs(coeff[4,0,0]))
    ql,qh=abs_range(coeff[5,0,0],total[5]-abs(coeff[5,0,0]))
    if pl*pl+ql*ql==0:raise ValueError("orientation needs subdivision")
    klo=2*sqrtlo(F(al*al,ph*ph+qh*qh))
    khi=2*sqrtup(F(ah*ah,pl*pl+ql*ql))
    full=klo>=F(1,80)
    # Cutoff: the diameter is 2/|k|. If not fully bounded use a truncated curve.
    if full:
        model=MODEL_LOWER*fourthlo(1/khi)
        tail=F(0)
    else:
        qlo=min(F(1),T*klo/2)
        z=qlo*qlo
        ser=sum(c*z**m for m,c in enumerate(TWEIGHTS))
        model=F(24)*(ser.numerator*2**40//ser.denominator)/2**40
        if khi*T>=2:model=min(model,MODEL_LOWER*fourthlo(1/khi))
        tail=F(40,81)*sqrtup(1-z)
    # Leave rational room for midpoint quadrature separately in the driver.
    target=model-tail
    g8=coeff[:8];bounds=total[:8];unit=max(bounds)
    if not unit:raise ValueError("zero shape")
    C=1<<SHAPE_BITS;gens=[];rest=bounds.copy();rnd=np.zeros(8,dtype=object)
    for k in range(order+1):
        for i in range(k+1):
            col=g8[:,i,k-i]
            if k and not np.any(col):continue
            rest-=abs(col)
            q=(col*C)//unit;rnd+=(q*unit!=col*C)
            gens.append([int(j) for j in q])
    rad=[int(a+b) for a,b in zip(-((-rest*C)//unit),rnd)]
    # Clip polynomial basis [0:4] and a constant. Use per-component intervals;
    # this filter alone never decides that the circle equation has real roots.
    divisor=den*H**(2*g)
    clip=[(F(1),F(1))]
    for idx,fac in [(8,-4),(9,4),(10,4)]:
        v=coeff[idx,0,0];err=total[idx]-abs(v)
        q=sorted([F(fac*(v-err),divisor),F(fac*(v+err),divisor)])
        clip.append(tuple(q))
    clip.append((F(-T*T),F(-T*T)))
    scale=1<<12
    cc=[];cr=[]
    for a,b in clip:
        cp=(a+b)*scale/2
        mid=cp.numerator//cp.denominator
        err=max(abs(a*scale-mid),abs(b*scale-mid))
        cc.append(mid);cr.append(-((-err.numerator)//err.denominator))
    return full,target,gens,rad,cc,cr,klo,khi


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("output")
    ap.add_argument("--nx",type=int,default=1024)
    ap.add_argument("--directions",type=int,default=8)
    ap.add_argument("--family",choices=("b","chord"),default="b")
    opts=ap.parse_args();build_source(opts.output,opts.nx,opts.directions,opts.family)

