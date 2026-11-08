"""Exact ancillary identities in the angular, differentiation, and curvature arguments."""
from fractions import Fraction as F
from math import factorial,comb
from pathlib import Path
import json

def shift(p,c):
    return [sum(p[j]*comb(j,k)*c**(j-k) for j in range(k,len(p))) for k in range(len(p))]
def deriv(p,j):
    if j>=len(p):return [F(0)]
    return [p[i]*F(factorial(i),factorial(i-j)) for i in range(j,len(p))]
def conv(a,b,d):
    return [sum(a[j]*b[i-j] for j in range(i+1) if j<len(a) and i-j<len(b)) for i in range(d+1)]
def audit():
    # Central differences represented with real shift w, after replacing
    # i*eta by w. The inverse becomes D/sinh(D). The substitution is formal.
    for degree in (3,9,19,32):
        poly=[F((-1)**i*(i*i+3),i+1) for i in range(degree+1)]
        w=F(1,23);a,b=shift(poly,w),shift(poly,-w)
        diff=[(i-j)/(2*w) for i,j in zip(a,b)]
        sinh=[F(0) if j%2 else F(1,factorial(j+1)) for j in range(degree+1)]
        inv=[F(1)]
        for j in range(1,degree+1):inv.append(-sum(sinh[i]*inv[j-i] for i in range(1,j+1)))
        assert conv(inv,sinh,degree)==[1]+[0]*degree
        q=[F(0)]*(degree+1)
        for j,c in enumerate(inv):
            for i,v in enumerate(deriv(diff,j)):q[i]+=w**j*c*v
        assert q[:-1]==deriv(poly,1) and q[-1]==0

    for m in (3,8,25):
        for k in (2,4,11):
            n=k+m;G=[F((-1)**j*(j+1),1+j*j) for j in range(m+1)]
            head=[G[j]*F(comb(n,j),comb(m,j)) for j in range(m+1)]
            for j in range(m+1):
                assert head[j]*F(factorial(n-j),factorial(m-j))*F(factorial(m),factorial(n))==G[j]
                assert head[j]*(n-j)/n==G[j]*F(comb(n-1,j),comb(m,j))

    # Quantitative tail: use n>=6. sqrt(n/pi)/n^2 <= 1/9 even at n>=3.
    assert 2+4*(F(2,3)+F(10,9))<10
    assert 5*(F(2,5)/(1-F(2,5))**2-F(2,5))<4
    assert 9*(F(10,3)+F(8,9))<39 and F(41,3)+4+39<64
    # Rephrase 2*3^(4/3)<9 and .5+8/128^(5/6)<.8 over Q.
    assert 8*3**4<9**3
    assert (F(8)/(F(4,5)-F(1,2)))**6<128**5
    assert F(1,128*256**2)>F(15,16)**256
    print("Exact finite-difference inversion, inverse-differentiation coefficients and signed-curvature bounds verified.")
    Path(__file__).with_name("n256_identity_audit.json").write_text(json.dumps(
       {"finite_difference_degrees":[3,9,19,32],"inverse_differentiation":"passed","signed_curvature_constants":"passed"},
       indent=2)+"\n")
if __name__=="__main__":audit()

