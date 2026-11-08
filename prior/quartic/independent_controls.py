"""Reconstruct the initial Bernstein tables by tensor interpolation on an integer grid."""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial, lcm
from pathlib import Path
import hashlib
import struct
import numpy as np


def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def mult(a,b):return a[0]*b[0]-a[1]*b[1],a[1]*b[0]+a[0]*b[1]
def factor(a,t):return a[0]*t,a[1]*t
def square(a):return mult(a,a)
def norm(a):return a[0]**2+a[1]**2
def conjugate(a):return a[0],-a[1]


def dirs8():
    # Four vectors and their perpendiculars, rotated by (264+23i)/265.
    base=[(F(1),F(0)),(F(12,13),F(5,13)),(F(119,169),F(120,169)),(F(5,13),F(12,13))]
    return [mult((F(264,265),F(23,265)),z) for z in base+[(-b,a) for a,b in base]]


def values(e,x,y):
    # Directly evaluate G_c(z)-G_c(-1) for c=0 and c=1, independent of factoring.
    z=factor(mult(e,(x,y)),F(5,2));z2=square(z);z3=mult(z2,z);z4=square(z2)
    f0=add(sub(z4,factor(z2,2)),(1,0))
    f1=sub(add(factor(z,4),(F(8,3),0)),factor(z3,F(4,3)))
    q0=factor(f0,F(2,5)**4);q1=factor(f1,F(2,5)**4)
    cross=mult(q1,conjugate(q0))
    return norm(q0),*cross,norm(q1),*q0,*q1,1


def polys():
    out=[]
    for j in range(9):
        p=[1]
        for k in range(9):
            if k==j:continue
            q=[0]*(len(p)+1)
            for i,v in enumerate(p):q[i]-=(k-4)*v;q[i+1]+=v
            p=q
        out.append([v*(-1)**(8-j)*comb(8,j) for v in p])
    return out


def make(path,n):
    # W_X and W_Y are scaled Lagrange interpolants. Use a blossom for BOTH
    # verifications of the Y controls (the primary program uses monomial nets).
    mat=polys();a=factorial(8)
    xs=np.array([[sum(c*p**k*n**(8-k) for k,c in enumerate(pol))
                  for pol in mat] for p in [4*(2*j+1-n) for j in range(n)]],dtype=object)
    m=lcm(*(comb(8,k) for k in range(9)))
    b=np.array([[sum(comb(j,l)*comb(8-j,k-l)*4**l*(-4)**(k-l)
                     for l in range(max(0,k-8+j),min(j,k)+1))*(m//comb(8,k))
                 for j in range(9)] for k in range(9)],dtype=object)
    y=np.array(mat,dtype=object)@b
    entries=dirs8()
    with Path(path).open("wb") as f:
        f.write(struct.pack("<4i",n,len(entries),36,4))
        for dr in entries:
            data=np.array([[values(dr,F(x-4),F(h-4)) for h in range(9)] for x in range(9)],dtype=object)
            result=[]
            for j in range(9):
                v=data[:,:,j];den=lcm(*(F(t).denominator for t in v.flat))
                integers=np.array([int(t*den) for t in v.flat],dtype=object).reshape(9,9)
                tensor=((xs@integers@y)*(1<<36))//(den*a*a*n**8*m)
                assert max(map(abs,tensor.flat))<2**60
                result.append(np.array(tensor,dtype="<i8"))
            np.stack(result,axis=1).tofile(f)
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


if __name__=="__main__":
    import sys
    print(make(sys.argv[1],int(sys.argv[2])))

