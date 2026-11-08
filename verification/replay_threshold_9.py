"""Replay all degree ranges for the sufficient threshold nine using standard Python."""
from pathlib import Path
import subprocess,sys,json
if __name__=='__main__':
    if not __debug__ or sys.flags.optimize:raise SystemExit('Assertions must be enabled.')
    g=Path(__file__).resolve().parent
    commands=[['replay_threshold_15.py'],['verify_conditional_line.py'],
              ['verify_conditional_line.py','--degrees','9','--L','41/50','--out','conditional_nonic_primary.json'],
              ['audit_conditional_line.py'],['verify_reciprocal_strips.py'],['joint_circle_certificate.py'],
              ['audit_low_degree_transfer.py']]
    for command in commands:subprocess.run([sys.executable,'-S','-B',*command],cwd=g,check=True)
    def data(name):return json.loads((g/name).read_text())
    assert sorted(map(int,data('conditional_independent.json')))==list(range(9,15))
    assert all(str(n) in data('joint_circle_intervals.json') for n in range(9,15))
    assert all(str(n) in data('low_degree_transfer_independent.json') for n in range(9,15))
    print('All analytic-range, finite-line, strip and chord verification programs for n>=9 passed.')

