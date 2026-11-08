"""Reconstruct every closed parameter interval of the quartic line check."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as F
from intervals import shape,BITS
import argparse,json

P=Path(__file__).resolve().parent
def build(N=128,bits=24,output_dir=None):
    target=P if output_dir is None else Path(output_dir)
    rows=[];meta=[]
    for k in range(3):
        for j in range(N):
            data,f,g,gamma=shape(F(j,N),F(j+1,N),k)
            nums=[]
            for z in data:
                for v in (z.re,z.im):nums.extend(v.bounds(bits))
            nums.extend((f.bounds(bits)[1],g.bounds(bits)[1]))
            ident=k*N+j
            rows.append(' '.join(map(str,[ident,*nums])))
            meta.append({'id':ident,'sheet':k,'left':str(F(j,N)),
                         'right':str(F(j+1,N)),'gamma_lower':str(gamma.lower())})
    inp=target/f'coefficients_{N}.txt';inp.write_text('\n'.join(rows)+'\n',encoding='ascii')
    (target/f'parameters_{N}.json').write_text(json.dumps({'intervals_per_sheet':N,'coefficient_bits':bits,'rows':meta},indent=1)+'\n',encoding='utf-8')
    print('Created',len(rows),'closed parameter intervals with outward coefficient bounds.')
    return inp
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--parts',type=int,default=128);o=a.parse_args();assert __debug__;build(o.parts)
