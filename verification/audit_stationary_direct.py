"""Independent direct polynomial-division check of the stationary coefficients.

The affine row formula is not used to construct A2. Instead A2 is
computed by differentiating the polynomial flow twice, and inverse germs
are constructed by their triangular defining equations over exact dual
rationals. These finite examples supplement, and do not prove, the
general identity established analytically in the paper.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q

class Dual:
    def __init__(self,a=0,b=0):self.a=Q(a);self.b=Q(b)
    def __add__(self,v):
        v=v if isinstance(v,Dual) else Dual(v)
        return Dual(self.a+v.a,self.b+v.b)
    __radd__=__add__
    def __neg__(self):return Dual(-self.a,-self.b)
    def __sub__(self,v):return self+-dual(v)
    def __rsub__(self,v):return dual(v)+-self
    def __mul__(self,v):
        v=dual(v);return Dual(self.a*v.a,self.a*v.b+self.b*v.a)
    __rmul__=__mul__
    def __truediv__(self,v):
        v=dual(v);assert v.a
        return Dual(self.a/v.a,(self.b*v.a-self.a*v.b)/v.a**2)
    def __rtruediv__(self,v):return dual(v)/self
    def __eq__(self,v):
        v=dual(v);return self.a==v.a and self.b==v.b

def dual(v):return v if isinstance(v,Dual) else Dual(v)

def add(a,b):
    c=[0]*max(len(a),len(b))
    for j,v in enumerate(a):c[j]+=v
    for j,v in enumerate(b):c[j]+=v
    return c

def mul(a,b,K=None):
    L=len(a)+len(b)-1
    if K is not None:L=min(L,K+1)
    c=[0]*L
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<L:c[i+j]+=x*y
    return c

def power(a,n,K):
    b=[Q(1)]
    for _ in range(n):b=mul(b,a,K)
    return b+[Q(0)]*(K+1-len(b))

def derivative(a):return [j*a[j] for j in range(1,len(a))]

def divide(a,b):
    a=list(map(Q,a));b=list(map(Q,b));out=[Q(0)]*max(1,len(a)-len(b)+1)
    for j in range(len(a)-len(b),-1,-1):
        out[j]=a[j+len(b)-1]/b[-1]
        for k,v in enumerate(b):a[j+k]-=out[j]*v
    assert all(v==0 for v in a[len(b)-1:])
    return out,a[:len(b)-1]

def inverse(p,K):
    n=len(p)-1;assert p[n]==1 and p[n-1]==0
    c=[Q(1)]
    for j in range(1,K+1):
        c.append(Q(0))
        value=power(c,n,j)[j]
        for i in range(n):
            shift=n-i
            if j>=shift:value+=p[i]*power(c,i,j-shift)[j-shift]
        c[j]=-value/n
    h=[Q(1)]
    for j in range(1,K+1):
        v=(1-j)*c[j]-sum(h[i]*h[j-i] for i in range(1,j))
        h.append(v/2)
    return c,h

def constant_after_inverse(a,c,h):
    total=0
    for i,v in enumerate(a):total+=v*mul(power(c,i,i),h,i)[i]
    return total

def verify():
    count=0
    for n in range(4,10):
        b=Q(n-1,2*n);sb=b*(b+1)/2;tb=b*b*(b-1)*(b+3)/4
        for seed in (1,2):
            q=[Q(((-1)**(j+seed))*(j+seed),n+2) for j in range(n-1)]+[Q(0),Q(1)]
            qdot,r=divide(mul(q,q),derivative(q))
            _,h= inverse([Dual(v,r[j] if j<len(r) else 0) for j,v in enumerate(q)],2*n)
            hn=h[n].a;eta=h[n].b;h2=h[2*n].a
            c=Q(seed,3);d=Q(1,2);delta=d**n
            p=[delta*q[j]/d**j for j in range(n+1)];p[0]+=c
            z,w=divide(mul(p,p),derivative(p))
            zdot,_=divide(add([2*v for v in mul(p,w)],[-v for v in mul(derivative(w),z)]),derivative(p))
            A1=[v/2 for v in derivative(z)]
            A2=add(add([v/4 for v in derivative(zdot)],[v/4 for v in mul(derivative(derivative(z)),z)]),
                   [v/8 for v in mul(derivative(z),derivative(z))])
            ci,hi=inverse(p,2*n)
            assert all(isinstance(v,Q) for v in hi)
            assert hi[n]==b*c+delta*hn
            assert constant_after_inverse(A1,ci,hi)==-b*hi[n]
            assert constant_after_inverse(A2,ci,hi)==tb*(c*c+2*delta*c*hn/b)-b*delta*delta*eta+sb*delta*delta*h2
            count+=1
    print('Direct rational polynomial-flow and inverse-germ checks:',count,'passed.')

if __name__=='__main__':verify()
