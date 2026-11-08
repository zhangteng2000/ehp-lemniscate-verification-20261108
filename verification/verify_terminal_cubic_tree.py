"""Build and replay one full exact or alternate terminal cubic chart from its index files."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from itertools import product
from fractions import Fraction as Q
from tempfile import TemporaryDirectory
import subprocess,sys,os,json,argparse
G=Path(__file__).resolve().parent
LEVELS=((8,16,32,256,128),(16,32,64,512,256),(32,64,128,1024,256),(64,128,256,1024,512))
def audit(n,k,scheme,workers):
    d=G/f'cubic_terminal_{n}';goals=json.loads((d/f'{k}_goals.json').read_text())
    assert set(goals)=={str(i) for i in range(32)};src=d/(scheme+'.cs');out=[]
    trig=G/('quintic_sine_dyadic.txt' if scheme=='primary' else 'quintic_sine_independent.txt')
    with TemporaryDirectory(dir=d) as tt:
        tmp=Path(tt);exe=tmp/'cubic.exe'
        if os.name=='nt':
            cc=os.environ.get('CSC') or str(Path(os.environ['WINDIR'])/'Microsoft.NET/Framework64/v4.0.30319/csc.exe')
            subprocess.run([cc,'/nologo','/checked+','/optimize+','/platform:x64','/out:'+str(exe),src],check=True)
        else:
            subprocess.run(['mcs','-checked+','-optimize+','-platform:x64','-out:'+str(exe),src],check=True)
        front=set(product(range(8),range(16),range(32)))
        for l,(q,a,b,T,H) in enumerate(LEVELS):
            if not front:break
            f=d/f'{k}{l}.txt';leaf=[tuple(map(int,x.split())) for x in f.read_text().splitlines()]
            assert set(leaf)<=front and len(set(leaf))==len(leaf)
            log=tmp/f'out{l}.txt'
            if leaf:subprocess.run([*([] if os.name=='nt' else ['mono']),exe,k,*map(str,(q,a,b,T,H)),f,log,trig,str(workers)],check=True)
            rows=[] if not leaf else [list(map(int,x.split())) for x in log.read_text().splitlines()]
            assert sorted(tuple(v[:3]) for v in rows)==sorted(leaf)
            for j,x,y,ct in rows:
                assert 0<=ct<=T*H
                need=max(Q(goals[str(i)]) for i in range(j*32//q,((j+1)*32+q-1)//q))
                assert Q((5 if k=='A' else 7)*333*ct,212*T*H)>=need,(n,k,l,j,x,y,ct)
            rest=front-set(leaf);out.append(len(leaf))
            front={(2*j+e,2*x+r,2*y+s) for j,x,y in rest for e,r,s in product(range(2),repeat=3)}
        assert not front
    print(n,k,scheme,'all terminal counts passed',out,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('n',type=int);p.add_argument('chart',choices=['A','B'])
    p.add_argument('scheme',choices=['primary','alternate']);p.add_argument('--workers',type=int,default=3)
    a=p.parse_args();assert __debug__;audit(a.n,a.chart,a.scheme,a.workers)

