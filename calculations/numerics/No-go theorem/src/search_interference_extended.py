"""Extended strain/power campaign, parallel numerical workers (not agents)."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from search_interference_temperature import optimize,OUT


def worker(args):
    d,o=args
    return optimize(d,o,20260922,zrange=(-4,5),maxiter=100)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    tasks=[(d,o) for d in (.5,1.683,5.,10.,30.,100.) for o in (.1,.3,1.,3.,10.)]
    rows=[]
    with ProcessPoolExecutor(max_workers=3) as pool:
        fs=[pool.submit(worker,a) for a in tasks]
        for f in as_completed(fs):
            r=f.result();rows.append(r)
            (OUT/'extended_search.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
            print(json.dumps({k:r[k] for k in ('strain_d_GHz','control_cap_GHz','witness_T_K',
                'local_success','certified_uniform_exclusion_above_K')}|
                {'certificate_passed':r['witness_certificate']['passed']}),flush=True)


if __name__=='__main__':main()
