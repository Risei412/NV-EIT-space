"""Gate 2.1: detuning-protocol audit and operational contours. See outputs/gate21_detuning_protocols.json, outputs/contour_protocols.json"""
# (executable version of the audited computation; identical to archived run)
import numpy as np,json,nv_model as nm
rows=[]
for (T,Bx) in [(300.,0.),(300.,1.),(110.,1.),(140.,1.),(80.,0.),(65.,0.)]:
    w,U=nm.dressed_ground((Bx,0,0.02)); H=nm.Hes((Bx,0,0.02))
    dp,dc=nm.legs(U); g=nm.gamma_oc_GHz(T,1.683)
    rows.append(dict(T=T,Bx=Bx,protocols={p:nm.scanC(H,dp,dc,g,p) for p in 'ABC'}))
import os
_OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                  '..','..','..','..','results','certificates')
json.dump(rows,open(os.path.join(_OUT,'gate21_detuning_protocols.json'),'w'),
          indent=1,default=float)
