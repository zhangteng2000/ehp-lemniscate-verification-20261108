"""Rebuild both complete quartic line covers and both circle checks.

Run from the code repository as
    python verification/quartic_line/verify_quartic_line.py --workers 2

Both coefficient inputs are regenerated on the three closed parameter
intervals; the fourth radical branch is the complex conjugate of the third.  Two
arithmetic circuits and two trigonometric enclosures are used; neither stored
counts nor a floating-point value can make an acceptance test pass.
"""

if not __debug__:
    raise RuntimeError('Verification requires assertions: do not run Python with -O or -OO.')

from pathlib import Path
from fractions import Fraction as F
import argparse,hashlib,json,os,shutil,subprocess,sys,tempfile

P=Path(__file__).resolve().parent
# During development the folder is next to github_workspace; in the published
# repository it is verification/quartic_line.
V=P.parent if (P.parent/'audit_quintic_rational.py').is_file() else P.parent/'github_workspace'/'verification'
sys.path.insert(0,str(V))

def compiler():
    choices=[os.environ.get('EHP_CSC'),os.environ.get('CSC'),shutil.which('csc'),shutil.which('mcs'),
             'C:/Windows/Microsoft.NET/Framework64/v4.0.30319/csc.exe']
    for q in choices:
        if q and Path(q).is_file():return q
    raise RuntimeError('Set EHP_CSC or CSC to csc or mcs supporting checked arithmetic.')

def check_cover(meta,counts,parts=128):
    data=json.loads(Path(meta).read_text(encoding='utf-8'))
    assert data['intervals_per_sheet']==parts and data['coefficient_bits']==24
    rows=data['rows'];assert len(rows)==3*parts
    ct={}
    for line in Path(counts).read_text(encoding='ascii').splitlines():
        i,c=map(int,line.split());assert i not in ct and 0<=c<=256*256
        ct[i]=c
    assert set(ct)==set(range(3*parts))
    minima=[]
    for sheet in range(3):
        value=[]
        for j in range(parts):
            row=rows[sheet*parts+j]
            assert row['id']==sheet*parts+j and row['sheet']==sheet
            assert F(row['left'])==F(j,parts) and F(row['right'])==F(j+1,parts)
            gamma=F(row['gamma_lower']);assert gamma>0
            bound=gamma*F(6*333*ct[row['id']],106*256*256)
            assert bound>F(29,20)
            value.append(bound)
        minima.append(min(value))
    return minima

def replay(workers=2):
    assert __debug__ and sys.flags.optimize==0
    from audit_quintic_rational import trig
    from build_inputs import build
    from independent_coefficients import build as independent_build
    from verify_circle_constants import verify as circle
    from independent_circle_audit import verify as independent_circle
    trig()
    csc=compiler()
    with tempfile.TemporaryDirectory(prefix='ehp-quartic-line-') as name:
        tmp=Path(name)
        source=tmp/'QuarticLineCubic.cs';shutil.copy2(P/source.name,source)
        exe=tmp/'quartic_line.exe'
        prefix='-' if Path(csc).stem.lower()=='mcs' else '/'
        flags=[prefix+q for q in ('nologo','checked+','optimize+','platform:x64','out:quartic_line.exe')]
        subprocess.run([csc,*flags,source.name],cwd=tmp,check=True)
        build(output_dir=tmp);independent_build(output_dir=tmp)
        result={}
        for scheme,inp,meta,sine in (
            ('expanded','coefficients_128.txt','parameters_128.json','quintic_sine_dyadic.txt'),
            ('factored','independent_coefficients_128.txt','independent_parameters_128.json','quintic_sine_independent.txt')):
            out=tmp/f'counts_{scheme}.txt'
            command=[str(exe),scheme,str(tmp/inp),str(out),'256','256',str(V/sine),str(workers)]
            if os.name!='nt':command=['mono',*command]
            subprocess.run(command,cwd=tmp,check=True)
            minima=check_cover(tmp/meta,out)
            rounded = ([F(61,25),F(37,25),F(73,50)] if scheme=='expanded'
                       else [F(49,20),F(38,25),F(151,100)])
            assert all(a>b for a,b in zip(minima,rounded))
            result[scheme]={'intervals':384,'sheet_minima':[str(q) for q in minima],
                            'minimum_gt':'29/20','coefficient_sha256':hashlib.sha256((tmp/inp).read_bytes()).hexdigest(),
                            'counts_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}
            shutil.copy2(out,P/f'recomputed_counts_{scheme}.txt')
        circle();independent_circle()
        result['source_sha256']={q.name:hashlib.sha256(q.read_bytes()).hexdigest()
            for q in sorted(P.iterdir()) if q.suffix in ('.py','.cs') and 'probe' not in q.name}
        result['status']='both complete closed line covers and both rational circle checks passed'
        (P/'quartic_line_audit.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(result['status'])
    return result

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--workers',type=int,default=2)
    o=a.parse_args();assert o.workers>=1;replay(o.workers)
