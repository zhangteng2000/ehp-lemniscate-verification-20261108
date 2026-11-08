"""Verify the listed rational shape implications without generating any proposed targets."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as Q
import json,sys,argparse
from verify_cubic_shape_implication import row
from verify_signed_source_moment import test
G=Path(__file__).resolve().parent
def run(n,kind):
    L={6:Q(73,100),7:Q(367,500)}[n];d=G/f'cubic_terminal_{n}';params=json.loads((G/'conditional_contours.json').read_text())
    gs=json.loads((d/f'{kind}_goals.json').read_text());sf=json.loads((d/f'{kind}_fallbacks.json').read_text())
    data=json.loads((G/f'proved_cubic_sufficient_profile_{n}.json').read_text())
    selected=[r for r in data if r[0]==kind];assert [r[1] for r in selected]==list(range(32))
    for rec in selected:
        k,i,z,goal=rec[:4];q=Q(gs[str(i)])
        assert Q(z)==Q(i+1,32)*(Q(5,16) if k=='A' else Q(9,4))
        if q>=Q(goal):
            ok,ans=row(n,L,k,Q(z),Q(goal),params)
            assert ok and ans==rec[4:],(n,k,i)
        else:
            f=G/sf[str(i)];s=json.loads(f.read_text());nn,kk,y,requested,mesh,closed,rr=s[:7]
            assert (nn,kk,y)==(n,k,i) and Q(requested)<=q and closed
            cnt=s[7] if len(s)>7 else 2
            ok,ans=test(n,L,k,i,Q(requested),tuple(mesh),f.name.startswith('paired_'),cnt,'_cond.' in f.name)
            assert ok and ans==rr,(n,k,i)
    print(n,kind,'all 32 scalar implications passed',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('n',type=int);p.add_argument('kind',choices=['A','B'])
    a=p.parse_args();assert __debug__;run(a.n,a.kind)

