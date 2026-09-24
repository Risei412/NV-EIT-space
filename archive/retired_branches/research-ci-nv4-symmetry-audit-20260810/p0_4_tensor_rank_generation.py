"""P0-4: tensor-rank generation from the four-orientation NV ensemble.

This follow-up corrects a subtlety in P0-3: at B=0 the ms=+1 and ms=-1
ground states are degenerate, so a one-sided spin-channel label is not a
basis-invariant observable under arbitrary transverse perturbations.

We therefore separate:
  (i) exact tetrahedral group moments/invariant dimensions,
  (ii) the spin-summed observable R_pm=(R_+ + R_-)/2,
  (iii) one-sided protocol-resolved cubic response.

Expected selection:
  rank 1: forbidden by tetrahedral averaging
  rank 2: only isotropic invariant
  rank 3: unique tetrahedral/octupolar invariant ~ Bx By Bz, accessible only
          in an appropriate protocol/parity sector
  rank 4: first even anisotropic invariant
"""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

import p0_3_nv4_orientation_symmetry_audit as p03

HERE = Path(__file__).resolve().parent
PRA = HERE.parent


def tetrahedral_axes():
    return np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float)/np.sqrt(3)


def invariant_dimension(n: int) -> int:
    """Molien coefficient for the proper tetrahedral rotation group T=A4."""
    idc = (n+1)*(n+2)//2
    c3 = 1 if n % 3 == 0 else 0
    c2 = sum(((-1)**k)*(k+1) for k in range(n+1))
    return int(round((idc + 8*c3 + 3*c2)/12))


def orbit_pm(B):
    _, t12 = p03.frame_sets()
    return 0.5*(p03.orbit_response(B, t12, "+1")
                + p03.orbit_response(B, t12, "-1"))


def hessian(func, h=1e-4):
    z = np.zeros(3); f0 = func(z); H = np.zeros((3,3))
    for i in range(3):
        e = np.zeros(3); e[i] = h
        H[i,i] = (func(e)+func(-e)-2*f0)/h**2
    for i in range(3):
        for j in range(i+1,3):
            ei=np.zeros(3); ej=np.zeros(3); ei[i]=h; ej[j]=h
            H[i,j]=H[j,i]=(func(ei+ej)-func(ei-ej)-func(-ei+ej)+func(-ei-ej))/(4*h*h)
    return f0,H


def rank4_estimate(r):
    ex = np.array([r,0,0.])
    nd = np.array([r,r,r])/np.sqrt(3)
    a = orbit_pm(ex); b = orbit_pm(nd)
    return dict(r_T=float(r), R_100=float(a), R_111=float(b),
                delta_R=float(a-b), c4=float((a-b)/((2/3)*r**4)))


def main():
    axes=tetrahedral_axes()
    M1=axes.sum(axis=0)
    M2=np.einsum('ai,aj->ij',axes,axes)
    M3=np.einsum('ai,aj,ak->ijk',axes,axes,axes)

    f0,H=hessian(orbit_pm,1e-4)
    eig=np.linalg.eigvalsh(H)

    parity=[]
    for B in (np.array([1e-3,1e-3,1e-3]),
              np.array([1e-3,-0.7e-3,0.4e-3]),
              np.array([2e-3,1e-3,-1.5e-3])):
        a=orbit_pm(B); b=orbit_pm(-B)
        parity.append(dict(B_T=B.tolist(), R_B=float(a), R_minusB=float(b),
                           relative_evenness_error=float(abs(a-b)/max(abs(a),abs(b),1e-300))))

    r4=[rank4_estimate(r) for r in (5e-4,1e-3,2e-3,3e-3)]
    rr=np.array([x['r_T'] for x in r4[:3]])
    cc=np.array([x['c4'] for x in r4[:3]])
    c4zero=float(np.polyfit(rr**2,cc,1)[1])

    _,t12=p03.frame_sets()
    cp=p03.dxyz(lambda B:p03.orbit_response(B,t12,"+1"),1e-4)
    cm=p03.dxyz(lambda B:p03.orbit_response(B,t12,"-1"),1e-4)

    out=dict(
        moments=dict(M1=M1.tolist(),M2=M2.tolist(),
                     M3_nonzero_value=float(4/(3*np.sqrt(3))),
                     cubic_contraction="sum_a(n_a.B)^3=(8/sqrt(3))BxByBz"),
        invariant_dimensions_rank_0_to_10=[invariant_dimension(n) for n in range(11)],
        spin_summed=dict(f0=float(f0),hessian=H.tolist(),eigenvalues=eig.tolist(),
                         parity=parity,rank4_samples=r4,small_field_c4=c4zero),
        one_sided_cubic=dict(dxyz_plus=float(cp),dxyz_minus=float(cm),
                             antisymmetric=float(0.5*(cp-cm)),
                             symmetric_residual=float(0.5*(cp+cm))),
        interpretation=("Tetrahedral averaging removes rank-1 and all traceless rank-2 response. "
                        "A protocol/parity sector may expose the unique rank-3 octupolar invariant; "
                        "for the basis-invariant spin-summed nearly field-even observable, rank-4 is "
                        "the first anisotropic magnetic sector."),
    )
    path=PRA/'results'/'tables'/'p0_4_tensor_rank_generation_reproduced.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
