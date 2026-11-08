"""Rebuild the all-degree EHP certificates. No floating design inputs or cached counters are accepted."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
import sys,subprocess,argparse,json
G=Path(__file__).resolve().parent
def call(args,site=False):
    subprocess.run([sys.executable,*([] if site else ['-S']),'-B',*map(str,args)],check=True,cwd=G)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=3)
    p.add_argument('--part',choices=['all','earlier','high','extension','archival-quartic'],default='all')
    a=p.parse_args();assert __debug__ and sys.flags.optimize==0
    if a.part=='all':
        call(['audit_cubic_literature.py'])
        call(['audit_common_identities.py'],True)
        call(['audit_inverse_coefficients.py'],True)
        call(['audit_stationary_identity.py'],True)
        call(['audit_stationary_direct.py'])
        call(['audit_capacity_majorants_exact.py'])
        call(['audit_octic_majorant_exact.py'])
        call(['audit_nonic_majorant_exact.py'])
        call(['uniform_circle_transfer.py'])
        call(['inverse_beta_moment_audit.py'])
    if a.part in ('all','earlier'):
        call(['quartic_line/verify_quartic_line.py','--workers',a.workers])
        call(['replay_low_packets.py','--degrees',5,'--workers',a.workers],True)
        call(['stationary_circle_transfer.py'])
        call(['independent_stationary_circle_transfer.py'])
    if a.part=='archival-quartic':
        call(['replay_low_packets.py','--degrees',4,'--workers',a.workers],True)
    if a.part in ('all','high'):
        call(['replay_threshold_8.py'])
        call(['single_coefficient_transfer.py'])
        call(['independent_single_coefficient_transfer.py'])
    if a.part in ('all','extension'):
        from audit_quintic_rational import trig
        trig()
        call(['audit_capacity_majorants_exact.py'])
        if a.part=='extension':
            call(['stationary_circle_transfer.py'])
            call(['independent_stationary_circle_transfer.py'])
        for n in (6,7):
            call([f'cubic_terminal_{n}/roundoff.py'])
            for k in 'AB':
                for s in ('primary','alternate'):call(['verify_terminal_cubic_tree.py',n,k,s,'--workers',a.workers])
                call(['verify_small_scalar_primary.py',n,k])
                call(['audit_small_scalar_profiles.py','--n',n,'--kind',k])
    print('All requested certificate ranges have been successfully recomputed.',flush=True)

