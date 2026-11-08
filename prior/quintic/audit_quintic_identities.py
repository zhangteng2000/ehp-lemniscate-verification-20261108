"""Polynomial-identity checks over Q by deterministic unisolvent grids."""
from fractions import Fraction as F
from itertools import combinations,product

def discriminant():
    tests=0
    for t in range(9):
        den=1+t*t;X=1-t*t;Y=2*t
        for h,b,c,r in product(range(-6,7),range(5),range(5),range(5)):
            # Clear den in the cubic coefficients and den^4 in the discriminant.
            a=Y;d1=3*h*X;d2=b*den-3*h*h*Y;d3=(c+h*r)*den-h**3*X
            disc=d1*d1*d2*d2-4*a*d2**3-4*d1**3*d3-27*a*a*d3*d3+18*a*d1*d2*d3
            d27=(27*(h*(Y*r-X*b)+Y*c)**2*den**2+4*Y*b**3*den**3
                 +108*h**3*X*c*den**3-36*h*h*b*b*den**4
                 +108*h**4*(X*r+Y*b)*den**3-108*h**6*den**4)
            assert disc+d27==0
            tests+=1
    # After the rational circle parametrization the degrees of the
    # cleared identity in t,h,b,c,r are at most 8,12,4,4,4.
    return tests

def shapes():
    for u,v in product(map(F,range(9)),repeat=2):
        pts=[F(0),F(1),u,v];mu=sum(pts)/4;z=[a-mu for a in pts]
        s=u+v-1;t=v-u
        assert -sum(a*a for a in z)/4==-(2+2*t*t+s*s)/16
        assert -sum(a**3 for a in z)/12==s*(1-t*t)/32
        def prim(a):return a**5-F(5,4)*(1+u+v)*a**4+F(5,3)*(u+v+u*v)*a**3-F(5,2)*u*v*a*a
        for i,j in combinations(range(4),2):
            a,b=pts[i],pts[j];r,l=[pts[k] for k in range(4) if k not in (i,j)]
            m=(a+b)/2;h=b-a
            assert prim(b)-prim(a)==-F(5,6)*h**3*((m-r)*(m-l)+h*h/20)
    # Each side for a fixed pair is a polynomial of individual degrees <=7.
if __name__=="__main__":
    assert __debug__
    n=discriminant();shapes()
    print(n,"discriminant evaluations on a determining grid, and all source identities passed.")

