#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q1 / DR-B -- chi_f scaling vs D2 at fixed criticality distance (RUNNER).
See README.md (preregistered). Pure NumPy; saves raw JSON incrementally."""
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    gammas=[0.5, 1.2, 1.4, 1.6, 1.8, 2.5],
    interior=[1.2, 1.4, 1.6, 1.8],
    sizes=[256, 512, 1024, 2048, 4096],
    n_real={256: 96, 512: 48, 1024: 24, 2048: 12, 4096: 6},
    window=0.25,
    fit_sizes=3,          # fit alpha and D2 over the largest `fit_sizes`
    seed=7,
    r1_tol=0.15,
    r2_min_delta=0.30,
)


def grp_matrix(N, gamma, rng):
    W = rng.standard_normal((N, N))
    W = (W + W.T) / np.sqrt(2.0)          # GOE: offdiag var 1, diag var 2
    H = N ** (-gamma / 2.0) * W
    H[np.diag_indices(N)] += rng.standard_normal(N)
    return H


def chi_window(E, U, win, V_diag, V_dense):
    """chi_f per window eigenstate for both perturbations."""
    n0, n1 = win
    cols = np.arange(n0, n1)
    dE = E[:, None] - E[None, cols]                        # (N, nwin)
    np.putmask(dE, np.abs(dE) < 1e-14, np.inf)             # m = n and exact degeneracies
    Mw_d = U.T @ (V_diag[:, None] * U[:, cols])            # (N, nwin)
    Mw_g = U.T @ (V_dense @ U[:, cols])
    chi_d = np.sum(Mw_d ** 2 / dE ** 2, axis=0)
    chi_g = np.sum(Mw_g ** 2 / dE ** 2, axis=0)
    return chi_d, chi_g


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": {k: (v if not isinstance(v, dict) else
                              {str(a): b for a, b in v.items()})
                          for k, v in cfg.items()}, "cells": {}}

    for gamma in cfg["gammas"]:
        for N in cfg["sizes"]:
            t0 = time.time()
            nr = cfg["n_real"][N]
            half = cfg["window"] / 2.0
            win = (int((0.5 - half) * N), int((0.5 + half) * N))
            ln_chi_d, ln_chi_g, ln_ipr, tail = [], [], [], []
            for _ in range(nr):
                H = grp_matrix(N, gamma, rng)
                E, U = np.linalg.eigh(H)
                ipr = np.sum(U[:, win[0]:win[1]] ** 4, axis=0)
                ln_ipr.extend(np.log(ipr).tolist())
                vd = rng.standard_normal(N)
                Vg = rng.standard_normal((N, N))
                Vg = (Vg + Vg.T) / (np.sqrt(2.0) * np.sqrt(N))
                cd, cg = chi_window(E, U, win, vd, Vg)
                ln_chi_d.extend(np.log(cd).tolist())
                ln_chi_g.extend(np.log(cg).tolist())
                tail.extend(cd.tolist())
            cell = dict(
                mln_chi_d=float(np.mean(ln_chi_d)),
                sln_chi_d=float(np.std(ln_chi_d) / np.sqrt(len(ln_chi_d))),
                mln_chi_g=float(np.mean(ln_chi_g)),
                mln_ipr=float(np.mean(ln_ipr)),
                n_states=len(ln_chi_d), n_real=nr,
                chi_d_q99=float(np.quantile(tail, 0.99)))
            results["cells"][f"{gamma}|{N}"] = cell
            print(f"[g={gamma:>3} N={N:>4}] ln chi_typ(diag) = "
                  f"{cell['mln_chi_d']:+.3f}, ln IPR = {cell['mln_ipr']:+.3f} "
                  f"({nr} real, {time.time()-t0:.1f} s)")
            with open(os.path.join(HERE, "results_raw.json"), "w") as f:
                json.dump(results, f, indent=1)

    # ---------------- fits ----------------
    fits = {}
    lnN_all = {N: np.log(N) for N in cfg["sizes"]}
    for gamma in cfg["gammas"]:
        Ns = cfg["sizes"][-cfg["fit_sizes"]:]
        x = np.array([lnN_all[N] for N in Ns])
        yd = np.array([results["cells"][f"{gamma}|{N}"]["mln_chi_d"] for N in Ns])
        yg = np.array([results["cells"][f"{gamma}|{N}"]["mln_chi_g"] for N in Ns])
        yi = np.array([results["cells"][f"{gamma}|{N}"]["mln_ipr"] for N in Ns])
        fits[str(gamma)] = dict(
            alpha_diag=float(np.polyfit(x, yd, 1)[0]),
            alpha_goe=float(np.polyfit(x, yg, 1)[0]),
            d2_meas=float(-np.polyfit(x, yi, 1)[0]),
            d2_theory=float(min(max(2.0 - gamma, 0.0), 1.0)))
    results["fits"] = fits

    # ---------------- preregistered verdicts ----------------
    inter = [fits[str(g)] for g in cfg["interior"]]
    dev = [abs(f["alpha_diag"] - (1.0 - f["d2_meas"])) for f in inter]
    n_ok = sum(d <= cfg["r1_tol"] for d in dev)
    alphas = [f["alpha_diag"] for f in inter]
    mono = all(np.diff(alphas) > 0) or all(np.diff(alphas) < 0)
    delta = abs(alphas[-1] - alphas[0])
    r1 = ("SUPPORTS" if n_ok >= 3 else
          ("PARTIAL (monotone tracking, offset > tol)"
           if mono and np.corrcoef(alphas, [1 - f["d2_meas"] for f in inter])[0, 1] > 0.9
           else "CHALLENGES"))
    r2 = ("SUPPORTS" if (mono and delta >= cfg["r2_min_delta"]) else
          "CHALLENGES")
    overall = ("SUPPORTS the D2 figure of merit"
               if (r1 == "SUPPORTS" and r2 == "SUPPORTS") else
               ("PARTIAL -- see R1/R2" if "CHALLENGES" not in (r1, r2)
                else "CHALLENGES the D2 figure of merit"))

    md = ["# Q1 / DR-B -- chi_f vs D2 (RESULTS)\n",
          f"gRP ensemble, sizes {cfg['sizes']}, window {cfg['window']:.0%}, "
          f"fits over the largest {cfg['fit_sizes']} sizes. Primary V = "
          f"diagonal; secondary V = GOE/sqrt(N).\n",
          "| gamma | alpha_diag | alpha_goe | D2_meas | D2_theory | "
          "1-D2_meas | dev |", "|---|---|---|---|---|---|---|"]
    for g in cfg["gammas"]:
        f_ = fits[str(g)]
        d = abs(f_["alpha_diag"] - (1 - f_["d2_meas"]))
        md.append(f"| {g} | {f_['alpha_diag']:+.3f} | {f_['alpha_goe']:+.3f} "
                  f"| {f_['d2_meas']:.3f} | {f_['d2_theory']:.2f} | "
                  f"{1-f_['d2_meas']:.3f} | {d:.3f} |")
    md.append(f"\n- **R1** (|alpha - (1-D2_meas)| <= {cfg['r1_tol']} at >= 3/4 "
              f"interior gammas): {n_ok}/4 within tol -> **{r1}**")
    md.append(f"- **R2** (monotone alpha, total change >= "
              f"{cfg['r2_min_delta']}): monotone = {mono}, delta = "
              f"{delta:.3f} -> **{r2}**")
    md.append(f"\n**Reading: {overall}.**")
    results["verdict"] = dict(R1=r1, R2=r2, overall=overall,
                              deviations=[float(d) for d in dev])
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
