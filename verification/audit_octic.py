"""Independent full-phase root-product certificate and octic circle comparison."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import argparse,json
from audit_peeled_line import PLO,PHI,root_floor,upper_sqrt,matrix_variance,atan_inv
from audit_conditional_line import check
from audit_low_degree_transfer import const,head
G=Path(__file__).resolve().parent
n=8;L=Q(39,50);t=Q(1,512);K=Q(53,100);S=1<<44
cap=[Q(j,1000) for j in (0,0,197,62,319,512,774,0,1000)]
def floor(y):return int((S*y).__floor__())
def ceil(y):return int((S*y).__ceil__())
def powup(x,d):
    d=Q(d)
    if d<0:return 1/root_floor(Q(x)**(-d.numerator),d.denominator)
    return root_floor(Q(x)**d.numerator,d.denominator)+Q(1,1<<76)
def jupper(v):
    # Split J(v)=int_0^(pi/2)sin(theta)^(v-1) dtheta at sin(theta)^2=1/2.
    s=z=Q(0);a=b=Q(1)
    for j in range(48):
        s+=a/(2**j*(2*j+v));z+=b/(2**j*(2*j+1))
        a*=Q(2*j+1,2*j+2);b*=Q(j+1-v/2,j+1)
    tail=Q(1,2**47)
    return (s+tail)*powup(2,-v/2)+(z+tail)*powup(2,Q(-1,2))
def nonlinear():
    stack=[((Q(0),Q(0),Q(0)),(Q(1),Q(1),Q(1)))];leaves=0
    while stack:
        a,b=stack.pop()
        if sum(a)>1:continue
        u,v,w=b
        M=(upper_sqrt((5*u+8*v)*(1-sum(a)))/4+3*w/8+
           5*upper_sqrt(3*w)*u/16+5*upper_sqrt(u)*v/8+11*u*u/128)
        if M<K:leaves+=1;continue
        j=max(range(3),key=lambda e:b[e]-a[e]);mid=(a[j]+b[j])/2;assert b[j]-a[j]>Q(1,4096)
        left=list(b);left[j]=mid;right=list(a);right[j]=mid
        stack.extend([(a,tuple(left)),(tuple(right),b)])
    return leaves
def curve():
    x,b,Il,Ih,gap,H=const(n);R=Q(3);c=Q(5,3);beta=Q(2,5)
    U=2*PHI*(1+(n-1)*Ih/4)
    e=t*R/root_floor(1-(t*R)**2,2)
    core=U*powup(R/2,x)*2*jupper(1-x)/PLO*(e+e*e)+(1+x)/(2+x)*t*t*powup(R,2+x)
    coef=(c/R)**4
    A=H*powup(R,x)*(c/R)**2/((2-x)*root_floor(1-coef,2))
    B=2*H**2*powup(c,4*x)*powup(R,-3*x)/(3*x*(1-coef))
    move=(4*A+2*B)*t/(beta-t)
    first=4*H*jupper(x)*powup(t,1-x)
    s=Q(2999,1000);assert s*s*(1+(t*R)**2)<R*R
    y=Q(ceil(powup(s,-x)),S);cut=208
    def kern(v,k):return 1+1/v+(1+PHI*k/2)/min(1,v-1)
    powers=[Q(1)]
    for i in range(2*cut+2):powers.append(Q(ceil(powers[-1]*y),S))
    # For an omitted ordered pair, degree h=i+j >= cut+2 and at most h pairs occur.
    h=cut+2;r=1-y
    tail=y**h*((2*h*h+(n+2)*h)/r+(4*h+n+2)*y/r**2+2*y*(1+y)/r**3)
    rows=[]
    for i in range(2,cut):
        rows.append(sum(kern(Q(i+j-1,n),abs(j-i)//n)*powers[i+j] for j in range(2,cut)
                        if (i-j)%n==0 and i+j>n+1))
    T=max(rows)+tail
    v=1/s;q=20
    a2=upper_sqrt(sum((kern(j-x,j)*v**j)**2 for j in range(2,q)))
    a2+=v**q*((3+2*q)/(1-v)+2*v/(1-v)**2)
    high=2*t*powup(R,1+x)*(H*H*T+2*H*a2)
    vals=[core,move,first,high];assert sum(vals)<Q(59,100)
    return [float(v) for v in vals]
def circle():
    x,b,Il,Ih,gap,H=const(n);entries=head(n,7);start=2*t/(1+t*t)
    stack=[(start,Q(3,4)),(Q(3,4),Q(1))];finished=0
    while stack:
        l,h=stack.pop();gain=Q(0)
        if h<=Q(3,4):
            a=root_floor(1-h*h,2);c=upper_sqrt(1-l*l)
            for m,(lg,hg,cs) in enumerate(entries,1):
                val=max(Q(0),lg*a**m-sum(s*h**j*c**(m-j) for j,s in cs.items()))**2/hg
                if m==1:val=max(val,(Il-1)*max(Q(0),b*a-K*h)**2)
                gain+=val
        if gap-2*PLO*gain<L*root_floor(l,n):finished+=1;continue
        assert h-l>Q(1,1<<29)
        p=(2*l+h)/3;q=(l+2*h)/3;stack.extend([(l,p),(p,q),(q,h)])
    return finished
def strips():
    m=n-1;nt=512;ng=128;nh=480;u=Q(13,25)
    @lru_cache(None)
    def sinlow(j,k):
        lo,hi=Q(j,nt),Q(k,nt)
        dist=0 if lo.__ceil__()<=hi else min(lo-lo.__floor__(),hi.__ceil__()-hi)
        a=PLO*dist/2;b=PHI*dist/2
        return floor(2*(a-b**3/6)*(1-b*b/2+a**4/24-b**6/720))
    # Horner evaluation of V,V',V'' with upper dyadic rounding after exact construction.
    data=[]
    for k in range(nh+1):
        h=Q(k,nh);v=w=z=Q(0)
        for j in range(m,-1,-1):
            z=z*h+2*w;w=w*h+v;v=v*h+(Q(j+1,n)*cap[j+1] if j else 0)
        data.append((ceil(v),ceil(w),ceil(z)))
    @lru_cache(None)
    def height(x,y):
        constant=(u*y/n).__floor__();lo,hi=0,nh;fac=m**m*x*x
        def passed(k):
            v,w,z=data[k];d=constant-v
            return d>0 and fac*d**(2*m-2)>(w*w+z*d)**m
        while hi-lo>1:
            k=(hi+lo)//2
            if passed(k):lo=k
            else:hi=k
        assert lo==0 or passed(lo)
        return lo
    angular=[sinlow(n*j,n*(j+1)) for j in range(nt)]
    areas=[]
    for b in range(ng):
        shift=b*(nt//ng);shift2=(b+1)*(nt//ng)
        areas.append(2*PLO*sum(height(x,sinlow(j+shift,j+1+shift2))
                               for j,x in enumerate(angular))/(nt*nh))
    assert min(areas)>Q(819,1000),float(min(areas))
    return min(areas)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--energy',default='octic_energy.json')
    p.add_argument('--contours',default='octic_contours.json');p.add_argument('--line-only',action='store_true');a=p.parse_args()
    assert __debug__;ll,hh=atan_inv(5);l0,h0=atan_inv(239);assert PLO<16*ll-4*h0<16*hh-4*l0<PHI
    c=check(n,L,json.loads((G/a.energy).read_text())['8'],json.loads((G/a.contours).read_text()))
    assert all(c[i]<=cap[i] for i in range(2,n+1))
    v=matrix_variance(n,cap,Q(13,25))/2;assert v<Q(97,100)
    if a.line_only:
        r={'variance_over_two':str(v),'strip_area':str(strips()),'status':'all line comparisons passed'}
        (G/'octic_line_independent.json').write_text(json.dumps(r,indent=2)+'\n')
        print('Independent octic line check: variance ratio',float(v),'; full-phase area',float(Q(r['strip_area'])))
        raise SystemExit(0)
    r={'variance_over_two':str(v),'strip_area':str(strips()),'nonlinear_leaves':nonlinear(),
       'short_arc_terms_display':curve(),'diameter_leaves':circle()}
    (G/'octic_independent.json').write_text(json.dumps(r,indent=2)+'\n')
    print('Independent variance ratio',float(v),'full phase area',float(Q(r['strip_area'])),
          'nonlinear leaves',r['nonlinear_leaves'],'curvature',sum(r['short_arc_terms_display']),
          'diameter leaves',r['diameter_leaves'])

