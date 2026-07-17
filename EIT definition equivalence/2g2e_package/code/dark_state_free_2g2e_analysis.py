"""
Reproduce the 2g+2e dark-state-free EIT analysis.

The model uses d_p=(1,0), d_c=(cos(theta),sin(theta)), so rank Omega=2
for sin(theta) != 0 and nonzero probe/control fields.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

gamma_e = 1.0
gamma_g = 0.02
Delta_e = 8.0
Omega_c = 0.8
theta = np.pi / 4

def response(delta, Gamma=0.0, path="none", theta_value=theta):
    delta = np.asarray(delta, dtype=complex)
    c, s = np.cos(theta_value), np.sin(theta_value)
    add1 = add2 = addg = 0.0
    if path == "full_optical":
        add1 = add2 = Gamma
    elif path == "all_coherences":
        add1 = add2 = addg = Gamma
    elif path == "selective_e2":
        add2 = Gamma
    elif path != "none":
        raise ValueError(path)
    a1 = gamma_e + add1 - 1j*delta
    a2 = gamma_e + add2 + 1j*Delta_e - 1j*delta
    gg = gamma_g + addg - 1j*delta
    beta = Omega_c**2/4
    den = gg*a1*a2 + beta*(c*c*a2 + s*s*a1)
    chi_full = (gg*a2 + beta*s*s)/den
    chi_cut = 1/a1
    return chi_full, chi_cut, chi_full-chi_cut

if __name__ == "__main__":
    d = np.linspace(-1.2, 1.2, 2401)
    full, cut, sector = response(d)
    print("Minimum Re chi_full:", full.real.min())
    print("rank[d_p,d_c] =", np.linalg.matrix_rank(
        np.array([[1.0, np.cos(theta)], [0.0, np.sin(theta)]])
    ))

    plt.figure()
    plt.plot(d, cut.real, label="cut")
    plt.plot(d, full.real, label="full")
    plt.xlabel("delta/gamma_e")
    plt.ylabel("normalized absorption")
    plt.legend()
    plt.tight_layout()
    plt.show()
