"""Exact scalar checks for the new quartic comparison. No parameter search."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as Q
import sys,json
P=Path(__file__).resolve().parent
V=P.parent if (P.parent/'audit_quintic_rational.py').is_file() else P.parent/'github_workspace'/'verification'
sys.path.insert(0,str(V))
from coherent_projection_transfer import constants
from diagonal_curvature_bound import certified
from audit_quintic_rational import PI

def verify():
    n=4;L=Q(29,20);t=Q(1,64);x=Q(1,4);b=Q(3,8)
    hn=Q(1,16)+3*t/8
    cc=constants(n)
    terms=certified(n,Q(5,2),t,Q(2,5),cc,hn_bound=hn,components=True)
    assert terms['total']<Q(7,5)<L
    assert (cc[4]-cc[2]).h<Q(139,200)
    delta=2*t/(1+t*t);mid=(1-t*t)/(1+t*t)
    assert delta==Q(128,4097) and mid==Q(4095,4097)
    h0=b*mid-delta/8
    assert h0==Q(12157,32776)
    assert Q(21,50)**4<delta
    margin=2*Q(333,106)*b*b*h0*h0+L*Q(21,50)-Q(139,200)
    assert PI.l>Q(333,106) and margin>Q(1,32)
    # f'(delta)>=L/4-2*pi*b^3>0 on (0,1/2], since -h'(delta)<1/2.
    assert b/Q(3,2)+Q(1,8)<Q(1,2)  # delta/sqrt(1-delta^2)<=1/sqrt(3)<2/3
    assert L/4>2*Q(22,7)*b**3
    # For delta>=1/2 the line bound alone gives a larger uniform gap.
    assert Q(6,5)**4>2
    assert (Q(29,24)-Q(139,200))/Q(6,5)>Q(1,32)
    # The short-arc loss gives a gap larger than 1/32 after rescaling.
    assert (L-Q(7,5))/Q(6,5)>Q(1,32)
    def up(v):return Q((v*(1<<40)).__ceil__(),1<<40)
    compact={k:up(v) for k,v in terms.items()}
    assert sum(compact[k] for k in ('core','motion','first','higher'))<Q(7,5)
    out={'curvature':{k:str(v) for k,v in compact.items()},'lambda':'29/20','tau':'1/64',
         'gap_upper':'139/200','middle_margin_lower':str(margin),
         'quartic_gap_constant':'1/32','status':'all rational acceptance tests passed'}
    (P/'circle_constants.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'curvature_components':{k:float(v) for k,v in terms.items() if k in ('core','motion','first','higher','total')},
                     'middle_margin_lower':str(margin),'gap_constant':'1/32'},indent=2))
if __name__=='__main__':assert __debug__;verify()
