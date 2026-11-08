"""Replay the extracted quartic and quintic code without manuscript files."""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from zipfile import ZipFile
from io import BytesIO
from hashlib import sha256
from tempfile import TemporaryDirectory
import sys,subprocess,argparse,os
G=Path(__file__).resolve().parent
def replay(n,outer,workers):
    name={4:'quartic_complete_verification_20261107.zip',5:'EHP_general_quintic_verified_20261107.zip'}[n]
    base=G.parent/'prior'/('quartic' if n==4 else 'quintic')
    if base.is_dir() and not outer:
        args=[sys.executable,'-B','verify_quartic.py' if n==4 else 'replay_general_quintic.py','--workers',str(workers)]
        if n==5:args.append('--second-source')
        subprocess.run(args,cwd=base,check=True)
        print('degree',n,'code-directory replay passed',flush=True)
        return
    if outer:
        with ZipFile(outer) as z:raw=z.read('prior/'+name)
    else:raw=(G.parent/'prior'/name).read_bytes()
    output=G/'low_replays';output.mkdir(exist_ok=True)
    with ZipFile(BytesIO(raw)) as packet,TemporaryDirectory(dir=output) as temp:
        assert packet.testzip() is None
        assert all(not Path(k).is_absolute() and '..' not in Path(k).parts for k in packet.namelist())
        packet.extractall(temp);base=Path(temp)/('quartic' if n==4 else 'certificate')
        cmd=[sys.executable,'-B','verify_quartic.py' if n==4 else 'replay_general_quintic.py','--workers',str(workers)]
        if n==5:cmd.append('--second-source')
        with (output/f'degree_{n}_output.txt').open('wb') as log:
            proc=subprocess.run(cmd,cwd=base,stdout=log,stderr=subprocess.STDOUT)
        proc.check_returncode()
    print('degree',n,'fresh packet replay passed',sha256(raw).hexdigest(),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--degrees',type=int,nargs='+',default=[4,5]);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--outer');a=p.parse_args();assert __debug__
    for n in a.degrees:replay(n,a.outer,a.workers)

