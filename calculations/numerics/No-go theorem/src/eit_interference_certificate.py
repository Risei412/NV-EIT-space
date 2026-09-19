"""Sufficient pole/residue EIT certificate for the passive 7-coherence model.

Not a necessary-and-sufficient definition for arbitrary multilevel systems.
The zero-detuning Schur complement equals nv_model.response exactly.
All probe-scanned coherences use the SAME -i*omega frequency convention.
"""
import numpy as np
import nv_model as nv

GEFF=2*6.3e-5+2e-6


def matrix(H,dc,z,gamma,Oc):
    M=np.zeros((7,7),complex)
    M[:6,:6]=gamma*np.eye(6)+2j*np.pi*(H-z*np.eye(6))
    M[:6,6]=1j*np.pi*Oc*dc
    M[6,:6]=1j*np.pi*Oc*dc.conj()
    M[6,6]=GEFF
    return M


def certificate(H,dp,dc,z,gamma,Oc,floor=.01):
    M=matrix(H,dc,z,gamma,Oc)
    f=np.r_[dp,0j]
    mu,V=np.linalg.eig(M)
    residues=(f.conj()@V)*(np.linalg.solve(V,f))
    slow=int(np.argmin(mu.real))
    a=float(mu[slow].real); center=float(mu[slow].imag)
    response=np.vdot(f,np.linalg.solve(M-1j*center*np.eye(7),f))
    G=np.linalg.inv(gamma*np.eye(6)+2j*np.pi*(H-z*np.eye(6))-1j*center*np.eye(6))
    Aoff=float(np.vdot(dp,G@dp).real)
    C=float((Aoff-response.real)/Aoff)
    reconstructed=np.sum(residues/(mu-1j*center))
    kappa=np.pi*Oc
    # Bauer-Fike disks of normal M(Oc=0): one isolated ground disk,
    # and six optical disks. Their real projections are disjoint here.
    separated=bool(gamma-GEFF>2*kappa)
    negative=bool(residues[slow].real < -1e-12)
    result=dict(passed=bool(separated and negative and C>=floor and Aoff>0),
        C_at_slow_center=C,Aoff=Aoff,Aon=float(response.real),
        slow_residue_real=float(residues[slow].real),
        slow_residue_imag=float(residues[slow].imag),
        slow_linewidth_rate=a,center_GHz=center/(2*np.pi),
        optical_min_decay_rate=float(np.min(np.delete(mu.real,slow))),
        decay_disk_margin=float(gamma-GEFF-2*kappa),
        no_ground_optical_pole_collision=separated,
        negative_narrow_residue=negative,
        pole_reconstruction_error=float(abs(reconstructed-response)),
        mode_condition_number=float(np.linalg.cond(V)))
    return result


def mixing_bound():
    # Strain and I_orb tensor H_g commute with ground-spin projectors.
    V0=nv.Hes((0,0,0),0)-np.kron(nv.I2,nv.Hgs((0,0,0)))
    return float(np.pi*np.ptp(np.linalg.eigvalsh(V0)))


def contrast_upper_bound(T,d,Oc):
    gamma=nv.gamma_oc_GHz(T,d)
    beta=(np.pi*Oc)**2
    return (mixing_bound()/gamma)**2*beta/(gamma*GEFF+beta)
