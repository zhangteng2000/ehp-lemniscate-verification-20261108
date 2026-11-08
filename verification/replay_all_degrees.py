"""Rebuild the all-degree EHP certificates. No floating design inputs or cached counters are accepted."""
from pathlib import Path
import sys,subprocess,argparse,json
G=Path(__file__).resolve().parent
def call(args,site=False):
    subprocess.run([sys.executable,*([] if site else ['-S']),'-B',*map(str,args)],check=True,cwd=G)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=3)
    p.add_argument('--part',choices=['all','earlier','high','extension'],default='all')
    a=p.parse_args();assert __debug__ and sys.flags.optimize==0
    if a.part in ('all','earlier'):call(['replay_low_packets.py','--workers',a.workers],True)
    if a.part in ('all','high'):call(['replay_threshold_8.py'])
    if a.part in ('all','extension'):
        from audit_quintic_rational import trig
        trig()
        for n,m,L in ((3,15,'6/5'),(6,20,'73/100'),(7,16,'367/500')):
            call(['verify_arc_capacity_transfer.py','--n',n,'--M',m,'--L',L])
        call(['audit_capacity_transfer.py'])
        for n in (6,7):
            call([f'cubic_terminal_{n}/roundoff.py'])
            for k in 'AB':
                for s in ('primary','alternate'):call(['verify_terminal_cubic_tree.py',n,k,s,'--workers',a.workers])
                call(['verify_small_scalar_primary.py',n,k])
                call(['audit_small_scalar_profiles.py','--n',n,'--kind',k])
    print('All requested certificate ranges have been successfully recomputed.',flush=True)

