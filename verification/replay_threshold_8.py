"""Reproduce the octic proof, with the earlier full threshold replay, from this directory."""
from pathlib import Path
import subprocess,sys
if __name__=='__main__':
    assert __debug__ and sys.flags.optimize==0
    wd=Path(__file__).resolve().parent
    for args in (['replay_threshold_9.py'],
                 ['verify_conditional_line.py','--degrees','8','--L','39/50','--parameters','octic_contours.json','--out','octic_energy.json'],
                 ['verify_root_product_strip.py'],['three_degree_chord_transfer.py','--n','8','--pure-jensen'],
                 ['audit_octic.py']):
        subprocess.run([sys.executable,'-S','-B',*args],cwd=wd,check=True)
    print('All constituent threshold checks, including both unrestricted octic checks, passed.')

