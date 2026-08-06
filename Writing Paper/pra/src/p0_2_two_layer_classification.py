"""P0-2: separate sector-induced signed response from spectroscopy labels.

This gate recomputes the full-Liouvillian T-B_perp grid and records two
independent classifications at every point:

1. sector_class: sign/detectability of the pathway-cut response
   Delta A = A_cut - A_full (mechanistic object),
2. spectral_class: four-model AIC/AICc classification of A_full(delta_2)
   among EIT, ATS, Fano and Lorentzian (spectroscopic object).

A point is called robust spectroscopic EIT only when

- EIT is the best of all four models,
- its IC advantage over the runner-up is at least ROBUST=6, and
- the legacy ATS-EIT difference also favors EIT by at least 6.

Thus a Fano-best spectrum can still carry positive sector-induced transparency,
but it is labelled ``Fano-shaped sector transparency``, never EIT.

Outputs
-------
Writing Paper/pra/results/tables/p0_2_two_layer_phase_diagram.csv
Writing Paper/pra/results/tables/p0_2_two_layer_summary.json
Writing Paper/pra/results/tables/p0_2_two_layer_phase_diagram.npz
Writing Paper/pra/results/figures/fig_p0_2_two_layer_phase_diagram.{png,pdf}
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PRA = HERE.parent
REPO = PRA.parents[1]
NOGO_SRC = REPO / "No-go theorem" / "src"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(NOGO_SRC))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap

import gate1_candidate_aic_bootstrap as g1  # noqa: E402
import p1_phase_diagram as p1               # noqa: E402
import run_prl_prediction as rp              # noqa: E402

ROBUST = float(g1.ROBUST)

SECTOR_CLASSES = (
    "sector_transparency",
    "sector_absorption",
    "sector_unresolved",
)
SPECTRAL_CLASSES = (
    "EIT",
    "ATS",
    "Fano",
    "Lorentzian",
    "spectral_ambiguous",
    "spectral_unresolved",
)
JOINT_CLASSES = (
    "spectroscopic_EIT",
    "Fano_shaped_sector_transparency",
    "Lorentzian_shaped_sector_transparency",
    "ATS_with_positive_sector_response",
    "sector_transparency_spectral_ambiguous",
    "control_induced_absorption",
    "unresolved",
)


def _bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _finite_float(value, default=float("nan")) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


def _ic_key(fit: dict) -> str:
    return "aicc" if bool(fit["used_aicc"]) else "aic"


def classify_spectrum(fit: dict) -> dict:
    """Four-model classification with an explicit runner-up margin."""
    key = _ic_key(fit)
    ics = {name: float(data[key]) for name, data in fit["models"].items()}
    ranked = sorted(ics, key=ics.get)
    best, runner = ranked[0], ranked[1]
    margin = float(ics[runner] - ics[best])
    d_ats_eit = float(fit["delta_aic_ats_eit"])

    robust_best = bool(margin >= ROBUST)
    if best == "EIT" and robust_best and d_ats_eit >= ROBUST:
        spectral = "EIT"
    elif best == "ATS" and robust_best and d_ats_eit <= -ROBUST:
        spectral = "ATS"
    elif best == "Fano" and robust_best:
        spectral = "Fano"
    elif best == "Lorentzian" and robust_best:
        spectral = "Lorentzian"
    else:
        spectral = "spectral_ambiguous"

    return {
        "spectral_class": spectral,
        "best_model": best,
        "runner_up_model": runner,
        "best_runner_margin": margin,
        "ic_used": key,
        "delta_ic_ats_eit": d_ats_eit,
        "legacy_eit_ats_verdict": str(fit["verdict"]),
        **{f"{name.lower()}_{key}": value for name, value in ics.items()},
        **{
            f"{name.lower()}_weight": float(fit["models"][name]["akaike_weight"])
            for name in fit["models"]
        },
    }


def classify_sector(cmax: float, dA: float, window_ok: bool, status: str) -> tuple[str, str]:
    """Mechanistic signed classification, independent of line-shape models."""
    if status != "ok" or not window_ok:
        return "sector_unresolved", "spectrum window unresolved"
    if not math.isfinite(cmax) or not math.isfinite(dA):
        return "sector_unresolved", "non-finite signed response"
    if abs(cmax) < p1.C_DETECT:
        return "sector_unresolved", f"|C|<{p1.C_DETECT:g}"
    if dA > 0:
        return "sector_transparency", "pathway sector lowers absorption"
    if dA < 0:
        return "sector_absorption", "pathway sector raises absorption"
    return "sector_unresolved", "zero signed response"


def joint_class(sector: str, spectral: str) -> str:
    if sector == "sector_absorption":
        return "control_induced_absorption"
    if sector != "sector_transparency":
        return "unresolved"
    return {
        "EIT": "spectroscopic_EIT",
        "ATS": "ATS_with_positive_sector_response",
        "Fano": "Fano_shaped_sector_transparency",
        "Lorentzian": "Lorentzian_shaped_sector_transparency",
        "spectral_ambiguous": "sector_transparency_spectral_ambiguous",
        "spectral_unresolved": "sector_transparency_spectral_ambiguous",
    }[spectral]


def evaluate_point(args) -> dict:
    T, Bx = map(float, args)
    t0 = time.time()
    row = {
        "T_K": T,
        "Bx_T": Bx,
        "Bz_T": float(p1.BZ0),
        "Oc_GHz": float(p1.OC),
    }
    try:
        d, Af, Ac, C, info = p1.adaptive_spectrum(T, Bx)
        status = str(info.get("status", "unknown"))
        row.update(status=status, half_MHz=float(info.get("half_MHz", np.nan)),
                   attempts=int(info.get("attempts", 0)))
        if status != "ok" or C is None or not np.all(np.isfinite(C)):
            row.update(
                Cmax=float("nan"), center_MHz=float("nan"), fwhm_MHz=float("nan"),
                Acut_at_peak=float("nan"), Afull_at_peak=float("nan"),
                dA_at_peak=float("nan"), at_edge=True, window_ok=False,
                sector_class="sector_unresolved",
                sector_reason=f"adaptive spectrum status={status}",
                spectral_class="spectral_unresolved",
                best_model="n/a", runner_up_model="n/a",
                best_runner_margin=float("nan"), ic_used="n/a",
                delta_ic_ats_eit=float("nan"), legacy_eit_ats_verdict="n/a",
                joint_class="unresolved",
            )
            return row

        cmax, center, fwhm, at_edge = p1._metrics(d, C)
        width = float(d[-1] - d[0])
        window_ok = bool((not at_edge) and np.isfinite(fwhm)
                         and fwhm <= p1.FILLS_WINDOW_FRAC * width)
        ipk = int(np.argmax(np.abs(C)))
        acut = float(Ac[ipk])
        afull = float(Af[ipk])
        dA = float(acut - afull)
        sector, sector_reason = classify_sector(cmax, dA, window_ok, status)

        fit = g1.fit_all(d, Af)
        spectral = classify_spectrum(fit)
        jc = joint_class(sector, spectral["spectral_class"])

        row.update(
            Cmax=float(cmax), center_MHz=float(center), fwhm_MHz=float(fwhm),
            Acut_at_peak=acut, Afull_at_peak=afull, dA_at_peak=dA,
            at_edge=bool(at_edge), window_ok=window_ok,
            sector_class=sector, sector_reason=sector_reason,
            **spectral, joint_class=jc,
        )
    except Exception as exc:  # pragma: no cover - preserved in output
        row.update(
            status=f"exception:{type(exc).__name__}",
            error=str(exc), Cmax=float("nan"), dA_at_peak=float("nan"),
            sector_class="sector_unresolved", spectral_class="spectral_unresolved",
            best_model="n/a", runner_up_model="n/a", best_runner_margin=float("nan"),
            joint_class="unresolved", window_ok=False, at_edge=True,
        )
    row["seconds"] = float(time.time() - t0)
    return row


def _write_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _upper_by_field(rows: list[dict], field_values: np.ndarray, klass_key: str,
                    klass_value: str) -> dict[str, float | None]:
    out = {}
    for B in field_values:
        vals = [r["T_K"] for r in rows
                if abs(r["Bx_T"] - float(B)) < 1e-10
                and r.get(klass_key) == klass_value]
        out[f"{B:.5f}"] = float(max(vals)) if vals else None
    return out


def _count_mixed_legacy(rows: list[dict], legacy_rows: dict[tuple[float, float], dict]) -> int:
    count = 0
    for row in rows:
        old = legacy_rows.get((row["T_K"], row["Bx_T"]))
        if old and old.get("klass") == "transparency" and old.get("best_model") != "EIT":
            count += 1
    return count


def load_legacy() -> dict[tuple[float, float], dict]:
    path = PRA / "results" / "tables" / "p1_phase_diagram.csv"
    with path.open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return {(_finite_float(r["T_K"]), _finite_float(r["Bx_T"])): r for r in rows}


def make_figure(rows: list[dict], Ts: np.ndarray, Bs: np.ndarray, out: Path) -> None:
    nt, nb = len(Ts), len(Bs)
    sector_order = list(SECTOR_CLASSES)
    spectral_order = list(SPECTRAL_CLASSES)
    joint_order = list(JOINT_CLASSES)
    sector = np.full((nt, nb), sector_order.index("sector_unresolved"), dtype=int)
    spectral = np.full((nt, nb), spectral_order.index("spectral_unresolved"), dtype=int)
    joint = np.full((nt, nb), joint_order.index("unresolved"), dtype=int)
    lookup = {(r["T_K"], r["Bx_T"]): r for r in rows}
    for i, T in enumerate(Ts):
        for j, B in enumerate(Bs):
            r = lookup[(float(T), float(B))]
            sector[i, j] = sector_order.index(r["sector_class"])
            spectral[i, j] = spectral_order.index(r["spectral_class"])
            joint[i, j] = joint_order.index(r["joint_class"])

    dB = float(Bs[1] - Bs[0]) if len(Bs) > 1 else 0.05
    dT = float(Ts[1] - Ts[0]) if len(Ts) > 1 else 5.0
    xe = np.append(Bs, Bs[-1] + dB) - dB / 2
    ye = np.append(Ts, Ts[-1] + dT) - dT / 2

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0), sharey=True)
    specs = [
        (axes[0], sector, sector_order,
         ["#2f6fbb", "#d65f5f", "#dddddd"],
         "mechanistic sector sign"),
        (axes[1], spectral, spectral_order,
         ["#2f6fbb", "#e69f00", "#7a5195", "#55a868", "#bbbbbb", "#eeeeee"],
         "four-model spectroscopy"),
        (axes[2], joint, joint_order,
         ["#2f6fbb", "#7a5195", "#55a868", "#e69f00", "#aaaaaa", "#d65f5f", "#eeeeee"],
         "joint label (no category mixing)"),
    ]
    for ax, data, labels, colors, title in specs:
        cmap = ListedColormap(colors)
        norm = BoundaryNorm(np.arange(-0.5, len(labels) + 0.5), cmap.N)
        ax.pcolormesh(xe, ye, data, cmap=cmap, norm=norm, shading="flat")
        ax.set_xlabel(r"$B_\perp$ (T)")
        ax.set_title(title, fontsize=9)
        handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="none") for c in colors]
        pretty = [s.replace("_", " ") for s in labels]
        ax.legend(handles, pretty, fontsize=5.5, loc="upper center",
                  bbox_to_anchor=(0.5, -0.19), ncol=2, frameon=False)
    axes[0].set_ylabel("temperature (K)")
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main(quick: bool = False, jobs: int = 2) -> dict:
    Ts, Bs = p1.grids(quick=quick)
    points = [(float(T), float(B)) for T in Ts for B in Bs]
    started = time.time()
    if jobs > 1:
        with Pool(jobs) as pool:
            rows = list(pool.imap(evaluate_point, points))
    else:
        rows = [evaluate_point(p) for p in points]

    tab = PRA / "results" / "tables"
    fig = PRA / "results" / "figures"
    tab.mkdir(parents=True, exist_ok=True)
    fig.mkdir(parents=True, exist_ok=True)
    _write_csv(tab / "p0_2_two_layer_phase_diagram.csv", rows)

    legacy = load_legacy()
    sector_counts = dict(Counter(r["sector_class"] for r in rows))
    spectral_counts = dict(Counter(r["spectral_class"] for r in rows))
    joint_counts = dict(Counter(r["joint_class"] for r in rows))
    mixed_new = [r for r in rows if r["joint_class"] == "spectroscopic_EIT"
                 and r["best_model"] != "EIT"]
    weak_eit = [r for r in rows if r["joint_class"] == "spectroscopic_EIT"
                and r["best_runner_margin"] < ROBUST]

    bx0_rows = [r for r in rows if abs(r["Bx_T"]) < 1e-12]
    bx0 = {
        "sector_counts": dict(Counter(r["sector_class"] for r in bx0_rows)),
        "spectral_counts": dict(Counter(r["spectral_class"] for r in bx0_rows)),
        "joint_counts": dict(Counter(r["joint_class"] for r in bx0_rows)),
        "has_positive_sector_response": any(
            r["sector_class"] == "sector_transparency" for r in bx0_rows),
        "has_robust_spectroscopic_EIT": any(
            r["joint_class"] == "spectroscopic_EIT" for r in bx0_rows),
    }

    candidate = min(rows, key=lambda r: abs(r["T_K"] - 70.0)
                    + 1000 * abs(r["Bx_T"] - float(rp.BX0)))
    gates = {
        "G_P02_1_two_layers_present": bool(sector_counts and spectral_counts),
        "G_P02_2_four_model_margin_enforced": not weak_eit,
        "G_P02_3_no_EIT_Fano_mixed_label": not mixed_new,
        "G_P02_4_Bperp_zero_audited": bool(bx0_rows),
        "G_P02_5_candidate_resolved": candidate["joint_class"] != "unresolved",
    }

    summary = {
        "what": "P0-2 two-layer full-Liouvillian phase classification",
        "definitions": {
            "sector_object": "Delta A = A_cut - A_full at the signed feature",
            "spectral_object": "A_full(delta_2) fitted to EIT/ATS/Fano/Lorentzian",
            "robust_spectral_rule": (
                "best of all four models and IC(best runner-up)>=6; "
                "EIT/ATS also require |IC_ATS-IC_EIT|>=6 with the matching sign"
            ),
        },
        "grid": {"T_K": [float(x) for x in Ts], "Bx_T": [float(x) for x in Bs]},
        "sector_counts": sector_counts,
        "spectral_counts": spectral_counts,
        "joint_counts": joint_counts,
        "sector_transparency_upper_T_by_field": _upper_by_field(
            rows, Bs, "sector_class", "sector_transparency"),
        "spectroscopic_EIT_upper_T_by_field": _upper_by_field(
            rows, Bs, "joint_class", "spectroscopic_EIT"),
        "fano_shaped_sector_transparency_upper_T_by_field": _upper_by_field(
            rows, Bs, "joint_class", "Fano_shaped_sector_transparency"),
        "legacy_mixed_transparency_Fano_count": _count_mixed_legacy(rows, legacy),
        "new_mixed_EIT_nonEIT_best_count": len(mixed_new),
        "Bperp_zero_audit": bx0,
        "candidate_70K": candidate,
        "gates": gates,
        "verdict": "PASS" if all(gates.values()) else "FAIL",
        "seconds": float(time.time() - started),
        "quick": bool(quick),
    }
    (tab / "p0_2_two_layer_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    # Machine-readable grids.
    lookup = {(r["T_K"], r["Bx_T"]): r for r in rows}
    sector_grid = np.array([[SECTOR_CLASSES.index(lookup[(float(T), float(B))]["sector_class"])
                             for B in Bs] for T in Ts], dtype=int)
    spectral_grid = np.array([[SPECTRAL_CLASSES.index(lookup[(float(T), float(B))]["spectral_class"])
                               for B in Bs] for T in Ts], dtype=int)
    joint_grid = np.array([[JOINT_CLASSES.index(lookup[(float(T), float(B))]["joint_class"])
                            for B in Bs] for T in Ts], dtype=int)
    np.savez(tab / "p0_2_two_layer_phase_diagram.npz",
             T_K=Ts, Bx_T=Bs, sector=sector_grid, spectral=spectral_grid,
             joint=joint_grid, sector_code=json.dumps(SECTOR_CLASSES),
             spectral_code=json.dumps(SPECTRAL_CLASSES),
             joint_code=json.dumps(JOINT_CLASSES))
    make_figure(rows, Ts, Bs, fig / "fig_p0_2_two_layer_phase_diagram")

    print(json.dumps({
        "sector_counts": sector_counts,
        "spectral_counts": spectral_counts,
        "joint_counts": joint_counts,
        "legacy_mixed_transparency_Fano_count": summary["legacy_mixed_transparency_Fano_count"],
        "Bperp_zero_audit": bx0,
        "candidate_joint_class": candidate["joint_class"],
        "gates": gates,
        "verdict": summary["verdict"],
    }, indent=2))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    main(quick=args.quick, jobs=args.jobs)
