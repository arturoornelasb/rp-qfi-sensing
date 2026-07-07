#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q1b -- enhancement-exponent gate, fresh seed (RUNNER).
See README.md (preregistered). Identical to Q1 except: seed = 8 and the
R1' gate compares (alpha_diag - alpha_goe) to (1 - D2_meas)."""
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
    fit_sizes=3,
    seed=8,
    r1_tol=0.15,
    r2_min_delta=0.30,
)


def grp_matrix(N, gamma, rng):
    W = rng.standard_normal((N, N))
    W = (W + W.T) / np.sqrt(2.0)
    H = N ** (-gamma / 2.0) * W
    H[np.diag_indices(N)] += rng.standard_normal(N)
    return H


def chi_window(E, U, win, V_diag, V_dense):
    n0, n1 = win
    cols = np.arange(n0, n1)
    dE = E[:, None] - E[None, cols]
    np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
    Mw_d = U.T @ (V_diag[:, None] * U[:, cols])
    Mw_g = U.T @ (V_dense @ U[:, cols])
    return (np.sum(Mw_d ** 2 / dE ** 2, axis=0),
            np.sum(Mw_g ** 2 / dE ** 2, axis=0))


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
            ln_chi_d, ln_chi_g, ln_ipr = [], [], []
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
            results["cells"][f"{gamma}|{N}"] = dict(
                mln_chi_d=float(np.mean(ln_chi_d)),
                mln_chi_g=float(np.mean(ln_chi_g)),
                mln_ipr=float(np.mean(ln_ipr)),
                n_states=len(ln_chi_d), n_real=nr)
            print(f"[g={gamma:>3} N={N:>4}] done ({time.time()-t0:.1f} s)")
            with open(os.path.join(HERE, "results_raw.json"), "w") as f:
                json.dump(results, f, indent=1)

    fits = {}
    for gamma in cfg["gammas"]:
        Ns = cfg["sizes"][-cfg["fit_sizes"]:]
        x = np.array([np.log(N) for N in Ns])
        yd = np.array([results["cells"][f"{gamma}|{N}"]["mln_chi_d"] for N in Ns])
        yg = np.array([results["cells"][f"{gamma}|{N}"]["mln_chi_g"] for N in Ns])
        yi = np.array([results["cells"][f"{gamma}|{N}"]["mln_ipr"] for N in Ns])
        fits[str(gamma)] = dict(
            alpha_diag=float(np.polyfit(x, yd, 1)[0]),
            alpha_goe=float(np.polyfit(x, yg, 1)[0]),
            d2_meas=float(-np.polyfit(x, yi, 1)[0]),
            d2_theory=float(min(max(2.0 - gamma, 0.0), 1.0)))
    results["fits"] = fits

    inter = [fits[str(g)] for g in cfg["interior"]]
    enh = [f["alpha_diag"] - f["alpha_goe"] for f in inter]
    dev = [abs(e - (1.0 - f["d2_meas"])) for e, f in zip(enh, inter)]
    n_ok = sum(d <= cfg["r1_tol"] for d in dev)
    alphas = [f["alpha_diag"] for f in inter]
    mono = all(np.diff(alphas) > 0) or all(np.diff(alphas) < 0)
    delta = abs(alphas[-1] - alphas[0])
    corr = float(np.corrcoef(enh, [1 - f["d2_meas"] for f in inter])[0, 1])
    r1 = ("SUPPORTS" if n_ok >= 3 else
          ("PARTIAL (monotone tracking, offset > tol)"
           if mono and corr > 0.9 else "CHALLENGES"))
    r2 = "SUPPORTS" if (mono and delta >= cfg["r2_min_delta"]) else "CHALLENGES"
    overall = ("SUPPORTS the D2 figure of merit (enhancement exponent)"
               if (r1 == "SUPPORTS" and r2 == "SUPPORTS") else
               ("PARTIAL -- see R1'/R2" if "CHALLENGES" not in (r1, r2)
                else "CHALLENGES the D2 figure of merit"))

    md = ["# Q1b -- enhancement-exponent gate, seed 8 (RESULTS)\n",
          "| gamma | alpha_diag | alpha_goe | enh = d-g | D2_meas | "
          "1-D2_meas | dev |", "|---|---|---|---|---|---|---|"]
    for g in cfg["gammas"]:
        f_ = fits[str(g)]
        e = f_["alpha_diag"] - f_["alpha_goe"]
        d = abs(e - (1 - f_["d2_meas"]))
        md.append(f"| {g} | {f_['alpha_diag']:+.3f} | {f_['alpha_goe']:+.3f} "
                  f"| {e:+.3f} | {f_['d2_meas']:.3f} | "
                  f"{1-f_['d2_meas']:.3f} | {d:.3f} |")
    md.append(f"\n- **R1'** (|enh - (1-D2_meas)| <= {cfg['r1_tol']} at >= 3/4 "
              f"interior): {n_ok}/4 within tol (corr {corr:.3f}) -> **{r1}**")
    md.append(f"- **R2**: monotone = {mono}, delta = {delta:.3f} -> **{r2}**")
    md.append(f"\n**Reading: {overall}.**")
    results["verdict"] = dict(R1p=r1, R2=r2, overall=overall,
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
