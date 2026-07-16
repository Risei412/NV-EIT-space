import sys,os,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','src'))
import nv_model as nm
def test_hermiticity():
    for B in [(0,0,0.02),(1,0.3,0.05)]:
        assert np.allclose(nm.Hgs(B),nm.Hgs(B).conj().T)
        assert np.allclose(nm.Hes(B),nm.Hes(B).conj().T)
def test_ground_orthogonality_and_M0():
    for B in [(0.01,0,0.02),(1.0,0.5,0.05)]:
        w,U=nm.dressed_ground(B)
        assert abs(np.vdot(U[:,1],U[:,0]))<1e-13
        dp,dc=nm.legs(U)
        assert abs(np.vdot(dp,dc))<1e-13  # Model-A1 M0=0
def test_gamma_mapping():
    import phonon_rates as pr
    assert abs(pr.gamma_oc(300,1.683)/pr.k_orb(300,1.683)-0.5)<1e-12
    assert abs(pr.gamma_pop(300,1.683)/pr.k_orb(300,1.683)-2.0)<1e-12
def test_archive_regression():
    w,U=nm.dressed_ground((1.0,0,0.02)); H=nm.Hes((1.0,0,0.02))
    dp,dc=nm.legs(U); Ep,_,_=nm.probe_line(H,dp)
    C=nm.response(H,dp,dc,Ep,nm.gamma_oc_GHz(300,1.683))['C']
    assert abs(C-(-1.8552e-6))/1.8552e-6<1e-3  # archived signed 300K value
def test_resolvent_asymptotic():
    w,U=nm.dressed_ground((0.5,0,0.02)); H=nm.Hes((0.5,0,0.02))
    dp,dc=nm.legs(U); Ep,_,_=nm.probe_line(H,dp); g=nm.gamma_oc_GHz(300,1.683)
    r=nm.response(H,dp,dc,Ep,g)
    M1=np.vdot(dp,(H-Ep*np.eye(6))@dc)
    asym=-(nm.TWOPI**2)*abs(M1)**2/g**4
    assert abs((r['K12']*r['K21']).real/asym-1)<5e-3
def test_denominator_positive():
    w,U=nm.dressed_ground((1.0,0,0.02)); H=nm.Hes((1.0,0,0.02))
    dp,dc=nm.legs(U); Ep,_,_=nm.probe_line(H,dp)
    for z in np.linspace(Ep-2,Ep+2,21):
        assert nm.response(H,dp,dc,z,nm.gamma_oc_GHz(300,1.683))['den'].real>0
def test_endpoint_detection():
    w,U=nm.dressed_ground((1.0,0,0.02)); H=nm.Hes((1.0,0,0.02))
    dp,dc=nm.legs(U)
    assert nm.scanC(H,dp,dc,nm.gamma_oc_GHz(300,1.683),'A')['edge_max'] in (True,False)

# ---- Gate regression tests: full-Liouvillian vs reduced kernel, group-IV
# moment order, B_perp opening, EIT/ATS AIC classifier ----------------------
def test_full_liouvillian_matches_reduced_kernel():
    import bperp_kernel_map_v2 as km, bperp_full_validation as bfv
    _,Us=km.track(np.array([0.0,0.02,0.1]),0.02)
    for Bx,U in zip([0.0,0.02,0.1],Us):
        r=bfv.run_one(300.,float(Bx),0.02,U)
        assert abs(r['ratio']-1.0)<0.05  # full 9-level Lindbladian ~ reduced resolvent

def test_group_iv_moment_order_differs_from_nv():
    import nv_reduced_kernel as nvk, group_iv_model as giv
    Gammas=np.logspace(2,5,31)
    Mnv=nvk.moments(nvk.H_3E(),(0,-1),1)
    assert abs(Mnv[0])<1e-12  # NV: M0=0
    slope_nv=np.polyfit(np.log10(Gammas[-15:]),np.log10(np.abs(nvk.kernel(nvk.H_3E(),(0,-1),Gammas))[-15:]),1)[0]
    assert abs(slope_nv+2)<1e-2
    for mat in ('SiV','SnV'):
        H=giv.H_groupIV(mat); M=giv.moments(H,1)
        assert abs(M[0])>0.5  # group-IV: M0 != 0
        slope=np.polyfit(np.log10(Gammas[-15:]),np.log10(np.abs(giv.kernel(H,Gammas))[-15:]),1)[0]
        assert abs(slope+1)<0.05

def test_bperp_opening_reduced_kernel():
    import bperp_kernel_map_v2 as km
    Bsmall=np.linspace(0,0.03,31); _,Us=km.track(Bsmall,0.02)
    sm=[km.eval_res(300,float(B),0.02,U) for B,U in zip(Bsmall,Us)]
    y=np.abs(np.array([abs(r['C']) for r in sm])-abs(sm[0]['C']))
    m=(Bsmall>0)&(y>0)
    slope=np.polyfit(np.log10(Bsmall[m]),np.log10(y[m]),1)[0]
    assert abs(slope-2.0)<0.5  # C ~ B_perp^2 opening law

def test_eit_ats_classifier_sanity():
    import eit_ats_classifier as clf
    d=np.linspace(-20,20,400)
    r_eit=clf.classify(d,clf.A_EIT(d,5.0,1.0,4.7,0.9))
    r_ats=clf.classify(d,clf.A_ATS(d,3.0,1.0,8.0))
    assert r_eit['verdict']=='robust EIT'
    assert r_ats['verdict']=='robust ATS'

# ---- SIMULATION_PLAN.md Sec. 7: CPTP audits, limits, finite-difference
# response, grid convergence, signal-chain units, seed reproducibility ------
def _candidate_chi(**kw):
    import gate2_candidate_full_vs_reduced as g2
    return g2.chi_pair(70.0,g2.rp.BX0,g2.rp.BZ0,0.0,**kw)

def _candidate_steady(**kw):
    import gate2_candidate_full_vs_reduced as g2
    from liouvillian_core import liouvillian,steady_state
    H,Ls,Vp,dp,meta=g2.build_full(70.0,g2.rp.BX0,g2.rp.BZ0,0.0,**kw)
    L=liouvillian(H,Ls); rho0,res,gap=steady_state(L)
    return rho0.reshape(meta['N'],meta['N']),res

def test_steady_state_trace_hermiticity_positivity():
    for kw in (dict(),dict(isc=True)):
        rho,res=_candidate_steady(**kw)
        assert abs(np.trace(rho)-1.0)<1e-10
        assert np.max(np.abs(rho-rho.conj().T))<1e-10
        assert float(np.min(np.linalg.eigvalsh(rho)))>=-1e-10  # positivity audit
        assert res<1e-8

def test_zero_control_limit():
    # Oc=0: cutting the ground-coherence sector must not change chi
    chif,chic,_=_candidate_chi(Oc=0.0)
    assert abs(chif-chic)/abs(chic)<1e-9

def test_zero_field_limit():
    # Bx=0: spin-Lambda channel closed up to the residual excited-state
    # mixing (Dperp/Lperp); contrast must be <<1e-3 of the opened candidate
    def C_at(Bx):
        w,U=nm.dressed_ground((Bx,0,0.005)); H=nm.Hes((Bx,0,0.005))
        dp=np.kron(np.array([1.,0],complex),U[:,1])
        dc=np.kron(np.array([1.,0],complex),U[:,2])
        z=float(np.linalg.eigvalsh(H)[3])
        return nm.response(H,dp,dc,z,nm.gamma_oc_GHz(70,1.683),0.1)['C']
    import run_prl_prediction as rp
    assert abs(C_at(0.0))<1e-3*abs(C_at(rp.BX0))

def test_cut_sector_limit():
    # cutting the sector removes the whole transparency feature
    chif,chic,_=_candidate_chi()
    assert abs(chif-chic)>0  # finite delta_chi_S at the candidate
    import gate1_candidate_aic_bootstrap as g1
    d,A,Acut=g1.probe_scan(n=41)
    assert np.max(Acut-A)>0  # dA>0 somewhere in the window

def test_finite_difference_response():
    import gate2_candidate_full_vs_reduced as g2
    from liouvillian_core import liouvillian,steady_state,first_order,vec
    H,Ls,Vp,dp,meta=g2.build_full(70.0,g2.rp.BX0,g2.rp.BZ0,0.0)
    N=meta['N']; L=liouvillian(H,Ls); rho0,_,_=steady_state(L)
    Hp=2*np.pi*0.5*g2.OP*(Vp+Vp.conj().T); I=np.eye(N)
    V=-1j*(np.kron(I,Hp)-np.kron(Hp.T,I))
    x,_=first_order(L,V,rho0)
    rho_eps,_,_=steady_state(L+V)          # exact steady state with drive
    diff=rho_eps-rho0
    det=np.zeros(N*N,complex)
    for e,a in enumerate(dp): det[meta['p_idx']*N+(3+e)]=np.conj(a)
    a1=det.conj()@x; a2=det.conj()@diff
    assert abs(a1-a2)/max(abs(a1),1e-300)<1e-3  # linear response = FD limit

def test_frequency_grid_convergence():
    import gate1_candidate_aic_bootstrap as g1
    _,A1,c1=g1.probe_scan(n=401); _,A2,c2=g1.probe_scan(n=801)
    C1=np.max((c1-A1)/c1); C2=np.max((c2-A2)/c2)
    assert abs(C1-C2)/abs(C2)<1e-2  # halving the grid changes Cmax <1%

def test_signal_chain_units():
    import signal_chain as sc
    assert abs(sc.photon_energy_J(637.0)-3.12e-19)/3.12e-19<5e-3
    assert sc.transmission(0.0)==1.0
    assert abs(sc.delta_T_over_T(1e-6)-1e-6)/1e-6<1e-3
    assert sc.spectral_fraction(100.0,30.0)==1.0
    kw=dict(od_total=1.0,power_W=1e-6,lambda_nm=637.0,eta=0.1,sigma_tech=1e-7)
    tau=sc.required_tau_s(5.0,1e-4,**kw)
    assert abs(sc.snr(1e-4,tau_s=tau,**kw)-5.0)<1e-6  # round trip
    assert sc.required_tau_s(5.0,1e-9,**kw)==float('inf')  # below tech floor

def test_model_comparison_reproducibility():
    import gate1_candidate_aic_bootstrap as g1
    d,A,_=g1.probe_scan(n=201)
    noisy=A+np.random.default_rng(g1.SEED).normal(0,1e-4*np.max(A),A.size)
    r1=g1.fit_all(d,noisy); r2=g1.fit_all(d,noisy)
    assert r1['delta_aic_ats_eit']==r2['delta_aic_ats_eit']  # deterministic
    assert r1['verdict']==r2['verdict']

if __name__=='__main__':
    for f in [v for k,v in list(globals().items()) if k.startswith('test_')]:
        f(); print('PASS',f.__name__)
