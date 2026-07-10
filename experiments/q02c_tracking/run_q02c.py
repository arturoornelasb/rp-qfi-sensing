#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q2c-lite -- does the sensitivity peak TRACK the crossover?
(preregistered in header, 2026-07-09). Fine W grid (2.0..3.5, step 0.25)
at L = 12 and L = 14 on the frozen Q2 machinery. W*(L) = argmax of a
quadratic fit to ln chi_typ(W); W_c(L) = interpolated <r~> = 0.46
crossing. FROZEN reading: TRACKING if the two drifts share sign and
|dW* - dW_c| <= 0.5; PINNED if |dW*| <= 0.25 while |dW_c| >= 0.5;
otherwise UNRESOLVED-AT-TWO-SIZES (L = 16 windowed extension stays
registered)."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "q02_platform_coexistence"))
from run_q02 import (CFG as Q2CFG, basis_half_filling, build_h,
                     dipole_diag, rt_from_levels)

HERE = os.path.dirname(os.path.abspath(__file__))
W_GRID = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5]
SIZES = {12: 200, 14: 64}


def main():
    t00 = time.time()
    cfg = Q2CFG
    rng = np.random.default_rng(77)
    results = {"cells": {}}
    curves = {}
    for L, n_real in SIZES.items():
        states, index = basis_half_filling(L)
        dim = len(states)
        vop = dipole_diag(L, states)
        half = cfg["window"] / 2
        w0, w1 = int((0.5 - half) * dim), int((0.5 + half) * dim)
        cols = np.arange(w0, w1)
        rs_c, chi_c = [], []
        for W in W_GRID:
            t0 = time.time()
            rts, ln_chi = [], []
            for _ in range(n_real):
                w = rng.uniform(-W, W, L)
                H = build_h(L, states, index, w, cfg["J"], cfg["V"])
                E, U = np.linalg.eigh(H)
                rts.extend(rt_from_levels(E[w0:w1]).tolist())
                dE = E[:, None] - E[None, cols]
                np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
                Mw = U.T @ (vop[:, None] * U[:, cols])
                ln_chi.extend(np.log(np.sum(Mw ** 2 / dE ** 2,
                                            axis=0)).tolist())
            r_m = float(np.mean(rts))
            c_m = float(np.mean(ln_chi))
            rs_c.append(r_m)
            chi_c.append(c_m)
            results["cells"][f"{L}|{W}"] = dict(r=r_m, mln_chi=c_m)
            print(f"[L={L} W={W}] r {r_m:.4f} lnchi {c_m:+.3f} "
                  f"({time.time()-t0:.0f} s)", flush=True)
        curves[L] = (np.array(rs_c), np.array(chi_c))

    md = ["# Q2c-lite -- W* tracking test (RESULTS)\n",
          f"W grid {W_GRID}; sizes {list(SIZES)}.\n",
          "| L | W* (chi peak) | W_c (r = 0.46) |", "|---|---|---|"]
    Ws = np.array(W_GRID)
    marks = {}
    for L in SIZES:
        rs, chis = curves[L]
        pc = np.polyfit(Ws, chis, 2)
        Wstar = float(np.clip(-pc[1] / (2 * pc[0]), Ws[0], Ws[-1])) \
            if pc[0] < 0 else float(Ws[np.argmax(chis)])
        Wc = float(np.interp(0.46, rs[::-1], Ws[::-1]))
        marks[L] = (Wstar, Wc)
        md.append(f"| {L} | {Wstar:.2f} | {Wc:.2f} |")
    dWs = marks[14][0] - marks[12][0]
    dWc = marks[14][1] - marks[12][1]
    md.append(f"\n- drifts 12 -> 14: dW* = {dWs:+.2f}, dW_c = {dWc:+.2f}")
    if abs(dWs) <= 0.25 and abs(dWc) >= 0.5:
        verdict = "PINNED: the sensitivity peak does not follow the drifting crossover."
    elif np.sign(dWs) == np.sign(dWc) and abs(dWs - dWc) <= 0.5:
        verdict = ("TRACKING: the sensitivity peak drifts WITH the "
                   "crossover -- the weaker-but-real coexistence holds; "
                   "the peak rides the transition's moving edge.")
    else:
        verdict = ("UNRESOLVED-AT-TWO-SIZES: drifts "
                   f"({dWs:+.2f} vs {dWc:+.2f}); the L = 16 windowed "
                   "extension stays registered.")
    md.append(f"\n**Reading: {verdict}**")
    results["marks"] = {str(L): marks[L] for L in SIZES}
    results["verdict"] = verdict
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
