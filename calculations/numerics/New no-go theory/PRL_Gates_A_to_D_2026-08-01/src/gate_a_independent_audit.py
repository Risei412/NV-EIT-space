"""Independent audit of NV-EIT PRL Gate A."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent

TWOPI = 2 * np.pi
s2 = 1 / np.sqrt(2)
Sz = np.diag([-1.0, 0.0, 1.0]).astype(complex)
Sx = s2 * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], complex)
Sy = s2 * np.array([[0, 1j, 0], [-1j, 0, 1j], [0, -1j, 0]], complex)
I3 = np.eye(3, dtype=complex)
sz_o = np.array([[-1, 0], [0, 1]], complex)
sx_o = np.array([[0, 1], [1, 0]], complex)
sy_o = np.array([[0, 1j], [-1j, 0]], complex)
I2 = np.eye(2, dtype=complex)

GE = 28.02495164
P = dict(Dpar=1.42, Lpar=5.33, Dperp=1.55 / 2, Lperp=0.20 / np.sqrt(2))

def hes(Bvec=(0.0, 0.0, 0.02), d=1.683, phi=0.0):
    Bx, By, Bz = Bvec
    dx, dy = d * np.cos(phi), d * np.sin(phi)
    H = np.kron(I2, P["Dpar"] * (Sz @ Sz - (2 / 3) * I3))
    H += -P["Lpar"] * np.kron(sy_o, Sz)
    H += P["Dperp"] * (np.kron(sz_o, Sy @ Sy - Sx @ Sx) - np.kron(sx_o, Sx @ Sy + Sy @ Sx))
    H += P["Lperp"] * (np.kron(sz_o, Sx @ Sz + Sz @ Sx) - np.kron(sx_o, Sy @ Sz + Sz @ Sy))
    H += dx * np.kron(sz_o, I3) - dy * np.kron(sx_o, I3)
    H += GE * (Bx * np.kron(I2, Sx) + By * np.kron(I2, Sy) + Bz * np.kron(I2, Sz))
    return H

H = hes()
dp = np.kron(np.array([0, 1], complex), np.array([0, 1, 0], complex))
dc = np.kron(np.array([1, 0], complex), np.array([1, 0, 0], complex))

def response(gamma, z=0.0, Oc=1.0, gg=6.3e-5):
    beta = (TWOPI * Oc) ** 2 / 4
    geff = 2 * gg + 2e-6
    G = np.linalg.inv(gamma * np.eye(6) + 1j * TWOPI * (H - z * np.eye(6)))
    K12 = np.vdot(dp, G @ dc)
    K21 = np.vdot(dc, G @ dp)
    S1 = np.vdot(dp, G @ dp)
    S2 = np.vdot(dc, G @ dc)
    dXi = -beta * K12 * K21 / (geff + beta * S2)
    Acut = float(np.imag(1j * S1))
    dA = float(-np.imag(1j * dXi))
    C = dA / Acut
    return K12, dXi, Acut, dA, C

def fit_nu(gamma, values, mask):
    x = np.log(gamma[mask])
    y = np.log(np.abs(values[mask]))
    return float(-np.polyfit(x, y, 1)[0])

gammas = np.logspace(2, 10, 300)
vals = [response(g) for g in gammas]
names = ["K12", "dXi", "Acut", "dA", "contrast"]
arrays = [np.array([v[i] for v in vals]) for i in range(5)]

tail = gammas > 1e8
physical_window = (gammas > 1e3) & (gammas < 2e4)
gamma_300K = 13349.549074059234
j300 = int(np.argmin(np.abs(gammas - gamma_300K)))

tail_exp = {n: fit_nu(gammas, a, tail) for n, a in zip(names, arrays)}
physical_exp = {n: fit_nu(gammas, a, physical_window) for n, a in zip(names, arrays)}
local_exp = {}
for n, a in zip(names[1:], arrays[1:]):
    local = -np.gradient(np.log(np.abs(a)), np.log(gammas))
    local_exp[n] = float(local[j300])

summary = {"gamma_300K": gamma_300K, "deep_tail_exponents": tail_exp, "physical_window_exponents": physical_exp, "local_exponents_near_300K": local_exp}

with open(OUT / "gate_a_independent_audit.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["Gamma", "abs_K12", "abs_dXi", "abs_Acut", "abs_dA", "abs_contrast"])
    for row in zip(gammas, *(np.abs(a) for a in arrays)):
        w.writerow(row)

with open(OUT / "gate_a_independent_audit.json", "w") as f:
    json.dump(summary, f, indent=2)

plt.figure(figsize=(7.2, 5.2))
plt.loglog(gammas, np.abs(arrays[3]), label="absolute signed absorption difference |dA|")
plt.loglog(gammas, np.abs(arrays[2]), label="baseline absorption |Acut|")
plt.loglog(gammas, np.abs(arrays[4]), label="normalized contrast |C|")
plt.axvline(gamma_300K, linestyle="--", label="Gamma(300 K)")
plt.xlabel("Gamma")
plt.ylabel("magnitude")
plt.title("Independent NV-EIT Gate A audit")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(OUT / "gate_a_independent_audit.png", dpi=180)
plt.close()

print(json.dumps(summary, indent=2))
