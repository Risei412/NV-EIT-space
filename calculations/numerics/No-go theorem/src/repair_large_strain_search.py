"""Independent starts for the large-strain low-power cases."""
import json
import numpy as np
from scipy.optimize import differential_evolution
import nv_model as nv
from search_interference_temperature import OUT,geometry
from refine_interference_campaign import refine


def main():
    results=[]
    for cap in (.1,.3):
        bounds=[(0,.5),(.005,.05),(0,np.pi/2),(0,np.pi/2),(-np.pi,np.pi),
                (-np.pi,np.pi),(0,2*np.pi),(-4,5),(0,cap),(0,1)]
        gamma=nv.gamma_oc_GHz(60.,100.)
        def f(x):
            H,dp,dc,z=geometry(x,100.)
            return -nv.response(H,dp,dc,z,gamma,x[8])['C']
        de=differential_evolution(f,bounds,seed=20260923,popsize=12,maxiter=180,
                 polish=False,integrality=[False]*9+[True])
        keys=['Bx','Bz','probe_theta','control_theta','probe_phase','control_phase',
              'strain_phi','z_fraction','Oc']
        start=dict(parameters=dict(zip(keys,map(float,de.x[:-1]))),channel=int(round(de.x[-1])))
        r=refine(100.,cap,start);results.append(r)
        (OUT/'large_strain_repair.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
        print(json.dumps(r),flush=True)


if __name__=='__main__':main()
