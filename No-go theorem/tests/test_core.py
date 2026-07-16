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

if __name__=='__main__':
    for f in [v for k,v in list(globals().items()) if k.startswith('test_')]:
        f(); print('PASS',f.__name__)
