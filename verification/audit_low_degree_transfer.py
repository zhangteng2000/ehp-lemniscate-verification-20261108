"""Independent rational exterior projections, circle coverings and unsymmetrized nonic strip test."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import json
from audit_peeled_line import PLO,PHI,root_floor,upper_sqrt,atan_inv
G=Path(__file__).resolve().parent
def const(n):
    x=Q(1,n);b=(1-x)/2;y=Q(1);Il=Q(1)
    for k in range(1,241):
        y*=Q(b+1-k,k);Il+=y*y
    Ih=Il+240*y*y
    root=root_floor(2,n);S=Q(0);T=Q(0);q=Q(1);v=Q(1)
    for k in range(66):
        S+=q/(2**k*(x+2*k));T+=v/(2**k*(2*k+1))
        q*=Q(2*k+1,2*k+2);v*=Q(1-x/2+k,k+1)
    C=2*root*(S/upper_sqrt(root+Q(1,1<<76))+T/upper_sqrt(2))
    return x,b,Il,Ih,2*PHI*(1+(n-1)*Ih/4)-C,upper_sqrt((n-1)*Ih/4)
@lru_cache(None)
def bino(t,k):
    r=Q(1)
    for z in range(k):r*=Q(t-z,z+1)
    return r
def head(n,count):
    x,b,Il,Ih,gap,H=const(n);out=[]
    for m in range(1,count+1):
        S=[];v=Q(0)
        for i in range(m):
            v+=(-1)**m*bino(b,i)*bino(-b,m-i);S.append(v)
        assert v+(-1)**m*bino(b,m)==0
        def inner(a):
            # Construct each Toeplitz correlation by its lag.
            d=sum(S[i]*a[i] for i in range(m));coef=Q(1)
            for lag in range(1,m):
                coef*=Q(lag-1-b,lag+b)
                d+=coef*sum(S[i]*a[i+lag]+S[i+lag]*a[i] for i in range(m-lag))
            vals=[d*Il-a[0]*S[0],d*Ih-a[0]*S[0]]
            return min(vals),max(vals)
        lo,hi=inner(S);assert lo>0;uu={}
        for j in range(1,m+1):
            for l in range(m-j+1):
                v=[Q(0)]*m;r=m-j-l
                for i in range(r+1):
                    h=(-1)**r*bino(b-j-l,i)*bino(-b-j,r-i)
                    for z in range(j):v[i+l+z]+=(-1)**z*comb(j-1,z)*h
                a,e=inner(v)
                w=H*Q((j+l)**(j+l),j**j*(l**l if l else 1))
                uu[j+l]=uu.get(j+l,0)+max(abs(a),abs(e))*w
        out.append((lo,hi,uu))
    return out
def circle(n):
    x,b,Il,Ih,gap,H=const(n);data=head(n,9 if n==9 else 4)
    L=Q(41,50) if n==9 else Q(7,8)
    stack=[(Q(1,128*n*n),Q(1,2)),(Q(1,2),Q(1))];count=0
    while stack:
        l,h=stack.pop();gain=Q(0)
        if h<=Q(1,2):
            cmin=root_floor(1-h*h,2);cmax=upper_sqrt(1-l*l)
            for m,(normlo,normhi,coef) in enumerate(data,1):
                err=sum(coef[k]*h**k*cmax**(m-k) for k in coef)
                gain+=max(Q(0),normlo*cmin**m-err)**2/normhi
            if n==9:gain=max(gain,(Il-1)*max(Q(0),b*cmin-Q(9,10)*h)**2)
        if gap-2*PLO*gain<L*root_floor(l,n):
            count+=1;continue
        assert h-l>Q(1,2**28),(n,float(l),float(h))
        m1=(2*l+h)/3;m2=(l+2*h)/3
        stack.extend([(l,m1),(m1,m2),(m2,h)])
    return count

def nonlinear():
    # Independent bounds using rational upper square roots.
    a, b = Q(11,40), Q(23,40)
    assert upper_sqrt(Q(10))/4 < 2*root_floor(a*b,2)
    assert upper_sqrt(Q(12))/8 < b
    assert max(5*upper_sqrt(Q(2))/24, 11*upper_sqrt(Q(2))/64) < b-a
    cross = Q(5,8)*upper_sqrt(Q(9,5))*upper_sqrt(Q(4,27))
    assert b+cross < Q(9,10)

def strips():
    m=8;n=9;U=Q(12,25)
    co=[Q(v,1000) for v in (0,0,103,14,84,311,521,785,0,1000)]
    A=[Q(j+1,n)*co[j+1] for j in range(m+1)]
    @lru_cache(None)
    def horner(h):
        # Discard constant and keep differentiated binomial coefficients.
        v=w=z=Q(0)
        for i in range(m,0,-1):
            z=h*z+w;w=h*w+v;v=h*v+A[i]
        return v*h,w*h+v,z*h+w
    def sinminor(a,b):
        # Linear distance from [a,b] to the integer lattice, followed by a cubic lower sine.
        j=a.__floor__();rr=a-j;ss=b-j
        d=Q(0) if ss>=1 else min(rr,1-ss)
        return max(Q(0),PLO*d-(PHI*d)**3/6)
    nt,ng,nh=336,36,384;minarea=Q(3)
    for g in range(ng):
        s=0
        for j in range(nt):
            lo,hi=Q(j,nt),Q(j+1,nt)
            xx=sinminor(n*lo,n*hi);yy=sinminor(lo+Q(g,ng),hi+Q(g+1,ng))
            if xx==0:continue
            def positive(b):
                h=Q(b,nh);v,w,z=horner(h);D=U*yy/n-v
                if D<=0:return False
                p=m*h*h*(m-xx*xx)+2*A[m-2]*xx
                return (m*xx*D)**2>p*(w*w+2*D*z)
            a,b=0,nh
            while b-a>1:
                jh=(a+b)//2
                if positive(jh):a=jh
                else:b=jh
            assert a==0 or positive(a);s+=a
        minarea=min(minarea,2*PLO*s/(nt*nh))
    assert minarea>Q(41,50),float(minarea)
    return float(minarea)
if __name__=='__main__':
    assert __debug__;l,h=atan_inv(5);ll,hh=atan_inv(239)
    assert PLO<16*l-4*hh<16*h-4*ll<PHI
    nonlinear();res={'nonic_independent_strip_area':strips(),'R9_upper':'9/10'}
    for n in range(9,15):res[str(n)]=circle(n);print(n,res[str(n)],flush=True)
    print('nonic independent strip area',res['nonic_independent_strip_area'])
    (G/'low_degree_transfer_independent.json').write_text(json.dumps(res,indent=1)+'\n', encoding='utf-8', newline='\n')

