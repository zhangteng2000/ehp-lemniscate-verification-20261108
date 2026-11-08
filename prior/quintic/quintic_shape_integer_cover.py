"""Integer-only exhaustive source-shape comparison for the quintic cubic-scale criterion."""
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
from hashlib import sha256
import argparse
import json
import numpy as np

BITS=23
S=1<<BITS
OFF=np.array(list(product(range(2),repeat=4)),dtype=np.int32)
PAIRS=list(combinations(range(4),2))


def I(x):return (int(x)*S,int(x)*S)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def neg(a):return -a[1],-a[0]
def times(a,k):return (a[0]*k,a[1]*k) if k>=0 else (a[1]*k,a[0]*k)
def div(a,d):return a[0]//d,-((-a[1])//d)
def mul(a,b):
    p=a[0]*b[0];q=a[1]*b[0];r=a[0]*b[1];s=a[1]*b[1]
    return np.minimum(np.minimum(p,q),np.minimum(r,s))//S,-((-np.maximum(np.maximum(p,q),np.maximum(r,s)))//S)
def sq(a):
    m=np.maximum(0,np.maximum(a[0],-a[1]));M=np.maximum(-a[0],a[1])
    return m*m//S,-((-M*M)//S)
def C(x):return I(x),I(0)
def plus(a,b):return add(a[0],b[0]),add(a[1],b[1])
def minus(a,b):return sub(a[0],b[0]),sub(a[1],b[1])
def cmul(a,b):return sub(mul(a[0],b[0]),mul(a[1],b[1])),add(mul(a[0],b[1]),mul(a[1],b[0]))
def square(a):return sub(sq(a[0]),sq(a[1])),times(mul(a[0],a[1]),2)
def rat(a,k,d):return div(times(a[0],k),d),div(times(a[1],k),d)
def abs_bounds(a):
    def raw(x):
        l=np.maximum(0,np.maximum(x[0],-x[1]));h=np.maximum(-x[0],x[1]);return l*l,h*h
    aa,bb=raw(a[0]),raw(a[1]);return aa[0]+bb[0],aa[1]+bb[1]
def det(a,b):return sub(mul(a[0],b[1]),mul(a[1],b[0]))


def classify(level,rows,ca=Fraction(7,4),cb=Fraction(5,2),unbounded=False):
    step=S>>level
    assert 0<step and rows.dtype==np.int32
    low=rows.astype(object if unbounded else np.int64)*step
    high=low+step
    u=((low[:,0],high[:,0]),(low[:,1],high[:,1]))
    v=((low[:,2],high[:,2]),(2*low[:,3]-S,2*high[:,3]-S))
    c=[C(0),C(1),u,v]
    # All variables are in [-1,1]; during the raw chord computation the
    # unconstrained distance rectangle has real and imaginary parts <=2.
    excluded=(low[:,0]>high[:,2])|(low[:,0]+low[:,2]>S)
    for i,j in PAIRS:
        excluded|=abs_bounds(minus(c[i],c[j]))[0]>S*S
    s=minus(plus(u,v),C(1));t=minus(v,u);ss=square(s);tt=square(t)
    f=rat(plus(plus(C(2),rat(tt,2,1)),ss),-1,16)
    g=rat(cmul(s,minus(C(1),tt)),1,32)
    dw=[]
    for i,j in PAIRS:
        a,b=c[i],c[j];k,l=[m for m in range(4) if m not in (i,j)]
        mid=rat(plus(a,b),1,2);d=minus(b,a);d2=square(d)
        w=cmul(cmul(d2,d),plus(cmul(minus(mid,c[k]),minus(mid,c[l])),rat(d2,1,20)))
        dw.append(rat(w,-5,6))
    for i,j in combinations(range(6),2):
        z=det(dw[i],dw[j]);excluded|=(z[0]>0)|(z[1]<0)
    d2=np.minimum(25*S*S//36+1,np.maximum.reduce([np.broadcast_to(abs_bounds(z)[1],len(rows)) for z in dw]))
    F2=np.asarray(abs_bounds(f)[0],dtype=object);G2=np.asarray(abs_bounds(g)[0],dtype=object)
    D2=np.asarray(d2,dtype=object)
    fpass=16*(16*ca.numerator)**20*F2**5>(17*ca.denominator)**20*D2**2*S**6
    gpass=64*(16*cb.numerator)**30*G2**5>(17*cb.denominator)**30*D2**3*S**4
    return np.where(excluded,1,np.where(fpass,2,np.where(gpass,3,0))).astype(np.uint8)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--depth",type=int,default=9)
    p.add_argument("--batch",type=int,default=8192)
    p.add_argument("--max-active",type=int,default=1600000)
    p.add_argument("--resume",type=Path)
    p.add_argument("--save",type=Path)
    p.add_argument("--log",type=Path)
    p.add_argument("--ca",default="7/4");p.add_argument("--cb",default="5/2")
    p.add_argument("--unbounded",action="store_true")
    args=p.parse_args()
    ca,cb=Fraction(args.ca),Fraction(args.cb)
    if args.resume:
        a=np.load(args.resume,allow_pickle=False);level=int(a["level"]);row=a["rows"]
        assert list(a["constants"])==[str(ca),str(cb)]
    else:level=0;row=np.zeros((1,4),dtype=np.int32)
    history=[]
    while True:
        counts=np.zeros(4,dtype=np.int64);nexts=[];dig=sha256()
        for start in range(0,len(row),args.batch):
            rr=row[start:start+args.batch];v=classify(level,rr,ca,cb,args.unbounded);dig.update(v.tobytes())
            counts+=np.bincount(v,minlength=4)
            nexts.append(rr[v==0])
        remaining=np.concatenate(nexts) if nexts else np.zeros((0,4),dtype=np.int32)
        record={"level":level,"counts":list(map(int,counts)),"classification_sha256":dig.hexdigest()}
        print(record,flush=True);history.append(record)
        if not len(remaining) or level>=args.depth or len(remaining)>args.max_active:
            break
        row=(2*remaining[:,None,:]+OFF[None,:,:]).reshape(-1,4);level+=1
    if args.save:np.savez_compressed(args.save,level=level,rows=remaining,constants=[str(ca),str(cb)])
    if args.log:args.log.write_text(json.dumps({"ca":str(ca),"cb":str(cb),"history":history,
              "residual_count":len(remaining)},indent=2)+"\n",encoding="utf-8")

if __name__=="__main__":main()
