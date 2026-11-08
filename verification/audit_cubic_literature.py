"""Exact audit of the elementary constants in Dahlke's cubic proof.

The source is DOI 10.5281/zenodo.20318611, Lemmas 4.1, 5.1--5.2
and Proposition 5.3. This is a source check, not a substitute for the
ray coverage and integral transformations in that paper.
"""
if not __debug__:
    raise RuntimeError('Verification requires assertions: do not use -O or -OO.')

from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

def arctan_bounds(inv,terms):
    s=sum(Q((-1)**j,(2*j+1)*inv**(2*j+1)) for j in range(terms))
    e=Q(1,(2*terms+1)*inv**(2*terms+1))
    return (s,s+e) if terms%2==0 else (s-e,s)

def polyadd(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,0)+v
    return {k:v for k,v in out.items() if v}

def polymul(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():out[i+j]=out.get(i+j,0)+u*v
    return {k:v for k,v in out.items() if v}

def scaled(a,c):return {k:v*c for k,v in a.items() if v*c}

def verify():
    a,b=arctan_bounds(5,36);c,d=arctan_bounds(239,10)
    pi_lo,pi_hi=Q(333,106),Q(355,113)
    assert pi_lo<16*a-4*d<16*b-4*c<pi_hi
    # The angle in the Machin identity is determined by tangent and range.
    t=Q(1,5);t2=2*t/(1-t*t);t4=2*t2/(1-t2*t2)
    assert (t4-Q(1,239))/(1+t4/Q(239))==1

    lo,hi=Q(14142135,10**7),Q(14142136,10**7)
    assert lo*lo<2<hi*hi
    lam_lo,lam_hi=lo-1,hi-1
    assert (1+lam_hi**2)**3<Q(1573,1000)**3*lam_lo
    assert 3*pi_hi/8<Q(1179,1000)
    sec_hi=Q(2613126,10**6)
    assert sec_hi**2>4+2*hi
    log_hi=Q(1615,1000)
    exponential_lower=sum(log_hi**j/factorial(j) for j in range(25))
    assert sec_hi+1+hi<exponential_lower
    # Cubing the Holder bound proves A(1/sqrt(2))<2.29 exactly.
    assert Q(1573,1000)**3*Q(1179,1000)*Q(1615,1000)**2<Q(229,100)**3
    assert pi_hi*pi_hi<Q(20,9)**2*2
    assert 3*Q(229,100)+Q(10,9)<8

    # General Laurent polynomial identities, with no sampled u or rho.
    x={1:Q(1,2),-1:Q(1,2)}
    T=polyadd(scaled(polymul(polymul(x,x),x),4),scaled(x,-3))
    assert T=={3:Q(1,2),-3:Q(1,2)}
    lhs=polymul(polyadd(T,{0:-1}),polymul({3:1,0:1},{3:1,0:1}))
    rhs=polymul(polyadd(T,{0:1}),polymul({3:1,0:-1},{3:1,0:-1}))
    assert lhs==rhs
    q={1:Q(1),0:Q(-1,2),-1:Q(1,2)}
    diff=polyadd(polymul(q,q),{2:Q(-1),1:Q(1),0:Q(-1)})
    assert diff=={0:Q(1,4),-1:Q(-1,2),-2:Q(1,4)}
    # In the final comparison, lambda^2<1 makes the remaining term
    # (-50 lambda^2+6 lambda^4)/7 strictly negative for lambda>0.
    assert 0<lam_lo<lam_hi<1 and -50+6*lam_hi**2<0
    receipt={'source':'Dahlke, DOI 10.5281/zenodo.20318611, pp. 4--8',
             'status':'exact constant and Laurent-polynomial checks passed',
             'lambda_enclosure':[str(lam_lo),str(lam_hi)],
             'A_upper':'229/100','B_upper':'10/9','K_upper':'8',
             'proof_by_degree_or_parameter_sampling':False,
             'scope':'The analytic branch coverage and integrals remain primary-text proof obligations; their source audit is separate.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('Cubic source audit: exact constants, Joukowski identity, and final strict sign passed.')

if __name__=='__main__':verify()
