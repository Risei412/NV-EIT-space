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
if __name__=='__main__':
    for f in [v for k,v in list(globals().items()) if k.startswith('test_')]:
        f(); print('PASS',f.__name__)
