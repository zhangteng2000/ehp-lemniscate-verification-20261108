"""Second exact evaluation of each terminal shape's coefficient and Crofton implications."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from math import comb
from audit_seeded_descent import passes
from audit_source_strips import source
from verify_peeled_line import variance
import argparse,json
G=Path(__file__).resolve().parent
def roundcaps(c,n):
    if c is None:return None
    for j in range(2,n-1):c[j]=Q((c[j]*10000).__ceil__(),10000)
    return c
def initial(n,kind,l,h,x,z):
    c=[None]*(n-1)+[Q(0),Q(1)];lower=[Q(0)]*(n+1)
    c[n-2]=Q(comb(n,2),3)*h*h*(1 if kind=='A' else z)
    c[n-3]=comb(n,3)*h**3*(z if kind=='A' else 1)
    lower[n-2]=Q(comb(n,2),3)*l*l*(1 if kind=='A' else x)
    lower[n-3]=comb(n,3)*l**3*(x if kind=='A' else 1)
    return c,lower
def one(n,kind,i,goal,record,alt):
    L={6:Q(73,100),7:Q(367,500)}[n];top={'A':Q(5,16),'B':Q(9,4)}[kind];x=i*top/32;z=(i+1)*top/32
    if Q(record[3])<=goal:
        q=Q(record[3]);c,_=initial(n,kind,0,L/q,x,z);c=roundcaps(passes(n,L,c,None),n);assert c
        u=Q(record[4]);assert variance(n,c,u)<2
        assert c[2:]==list(map(Q,record[5]))
        area=source(n,1,c,u,0,c[n-2],'unsigned',(672,1,18,768,690))
        assert area>L,(n,kind,i,'unsigned',float(area))
        return ['unsigned',float(area)]
    f=G/alt;v=json.loads(f.read_text());nn,k,j,q,mes,ok,assertions=v[:7]
    assert (nn,k,j)==(n,kind,i) and ok and Q(q)<=goal
    parts=v[7] if len(v)>7 else 2;cond='_cond.' in f.name;paired=f.name.startswith('paired_');typ='paired' if paired else 'signed'
    levels=[Q(j,100) for j in ((0,80,88,92,95,97,99,100) if cond else (0,80,90,95,100))]
    rows={(Q(r[0]),Q(r[1]),r[2]):r for r in assertions};answers=[]
    # A separate finer mesh and a sine table reconstructed from dyadic half angles.
    K=max(840,mes[0]);A=max(10,mes[1]);P=max(144,mes[2]);N=max(560,mes[3]);Hmax=(N*mes[4]+mes[3]-1)//mes[3]
    for a,b in zip(levels,levels[1:]):
        c,lower=initial(n,kind,a*L/Q(q),b*L/Q(q),x,z)
        c=roundcaps(passes(n,L,c,lower if cond else None),n)
        if c is None or lower[n-2]>c[n-2]:continue
        for part in range(parts):
            cc=c.copy();bot=cc[2]*part/parts;cc[2]*=Q(part+1,parts)
            if cond:
                low=lower.copy();low[2]=bot;cc=passes(n,L,cc,low,1)
                if cc is None or any(cc[j]<low[j] for j in range(2,n-1)):continue
            rr=rows[a,b,part];u=Q(rr[3]);a1=Q(0)
            if u and variance(n,cc,u)<2:a1=source(n,1,cc,u,lower[n-2],c[n-2],typ,(K,A,P,N,Hmax))
            if a1<=L:a1=max(a1,source(n,2,cc,bot,lower[n-2],c[n-2],typ,(K,A,P,N,Hmax)))
            assert a1>L,(n,kind,i,a,b,part,float(a1))
            answers.append(float(a1))
    return [typ,min(answers,default=1.0)]
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,required=True);p.add_argument('--kind',required=True)
    p.add_argument('--bins',nargs='+',type=int);a=p.parse_args();n=a.n;kind=a.kind;td=G/f'cubic_terminal_{n}'
    goals=json.loads((td/f'{kind}_goals.json').read_text());orig=json.loads((G/f'proved_cubic_sufficient_profile_{n}.json').read_text())
    special=json.loads((td/f'{kind}_fallbacks.json').read_text());out=[]
    for i in a.bins if a.bins else range(32):
        row=next(r for r in orig if r[0]==kind and r[1]==i)
        v=one(n,kind,i,Q(goals[str(i)]),row,special.get(str(i)));out.append((i,v))
        print(n,kind,i,*v,flush=True)
    tail='' if a.bins is None else '_'+str(out[0][0])+'_'+str(out[-1][0])
    (td/f'{kind}_scalar_audit{tail}.json').write_text(json.dumps(out)+'\n')

