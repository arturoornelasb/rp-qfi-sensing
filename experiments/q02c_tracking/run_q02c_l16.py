#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q2c-full -- the L = 16 third point of the W* tracking test
(preregistered in header, 2026-07-12; the registered "windowed-chi
extension" of Q2c-lite, and the three-size follow-up the published
paper names). Same frozen Q2 machinery and W grid (2.0..3.5, 0.25);
L = 16 half-filling (dim 12,870) with 6 disorder realizations per W
(window ~3.2k states/realization -> ~19k chi samples per W, ample).
INSTRUMENT NOTE (mandatory at this scale): numpy's dsyevd eigh is
BROKEN near N ~ 1e4 (repo pitfall, Q3); ALL solves use scipy
driver='evr'; realization #1 of every W cell is certified against an
independent eigvalsh reference (agreement <= 1e-12 x width) -- the Q2
policy. FROZEN reading over the three sizes (12/14 loaded from
Q2c-lite's frozen results): TRACKING-CONFIRMED if both size steps have
sign(dW*) = sign(dW_c) and cumulative |dW* - dW_c| <= 0.5; PINNED if
|W*(16) - W*(12)| <= 0.25 while |W_c(16) - W_c(12)| >= 0.75;
UNRESOLVED otherwise."""
import json
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh as seigh

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "q02_platform_coexistence"))
from run_q02 import (CFG as Q2CFG, basis_half_filling, build_h,
                     dipole_diag, rt_from_levels)

HERE = os.path.dirname(os.path.abspath(__file__))
W_GRID = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5]
L, N_REAL = 16, 6
CERT_TOL = 1e-12


def wstar_wc(ws, lnchi, rs):
    a = np.polyfit(ws, lnchi, 2)
    w_pk = float(-a[1] / (2 * a[0])) if a[0] < 0 else float(
        ws[int(np.argmax(lnchi))])
    w_pk = min(max(w_pk, ws[0]), ws[-1])
    wc = None
    for i in range(len(ws) - 1):
        if (rs[i] - 0.46) * (rs[i + 1] - 0.46) <= 0:
            t = (0.46 - rs[i]) / (rs[i + 1] - rs[i])
            wc = float(ws[i] + t * (ws[i + 1] - ws[i]))
            break
    return w_pk, wc


def main():
    t00 = time.time()
    cfg = Q2CFG
    rng = np.random.default_rng(78)
    states, index = basis_half_filling(L)
    dim = len(states)
    vop = dipole_diag(L, states)
    half = cfg["window"] / 2
    w0, w1 = int((0.5 - half) * dim), int((0.5 + half) * dim)
    cols = np.arange(w0, w1)
    print(f"[setup] L={L} dim={dim}, window {w1-w0}", flush=True)

    results = {"L16_cells": {}}
    rs_c, chi_c = [], []
    for W in W_GRID:
        rts, ln_chi = [], []
        for k in range(N_REAL):
            t0 = time.time()
            w = rng.uniform(-W, W, L)
            H = build_h(L, states, index, w, cfg["J"], cfg["V"])
            E, U = seigh(H, driver="evr")
            if k == 0:
                Eref = np.linalg.eigvalsh(H)
                width = float(Eref[-1] - Eref[0])
                agree = float(np.max(np.abs(E - Eref)))
                ok = agree <= CERT_TOL * width
                results["L16_cells"][f"cert|{W}"] = dict(
                    agree=agree, width=width, ok=bool(ok))
                print(f"  [cert W={W}] agree {agree:.2e} "
                      f"({'OK' if ok else 'FAIL'})", flush=True)
                if not ok:
                    raise SystemExit(f"INSTRUMENT-FAIL at W={W}")
            rts.extend(rt_from_levels(E[w0:w1]).tolist())
            dE = E[:, None] - E[None, cols]
            np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
            Mw = U.T @ (vop[:, None] * U[:, cols])
            ln_chi.extend(np.log(np.sum(Mw ** 2 / dE ** 2,
                                        axis=0)).tolist())
            del H, U, dE, Mw
            print(f"  [W={W} r{k+1}/{N_REAL}] {time.time()-t0:.0f}s "
                  f"(total {time.time()-t00:.0f}s)", flush=True)
        r_m, c_m = float(np.mean(rts)), float(np.mean(ln_chi))
        rs_c.append(r_m)
        chi_c.append(c_m)
        results["L16_cells"][f"{L}|{W}"] = dict(r=r_m, mln_chi=c_m)
        with open(os.path.join(HERE, "results_l16.json"), "w") as f:
            json.dump(results, f, indent=1)

    ws = np.array(W_GRID)
    w16, wc16 = wstar_wc(ws, np.array(chi_c), np.array(rs_c))
    prev = json.load(open(os.path.join(HERE, "results_raw.json")))
    w12, wc12 = prev["marks"]["12"]
    w14, wc14 = prev["marks"]["14"]
    results["summary"] = {"12": dict(W_star=w12, W_c=wc12),
                          "14": dict(W_star=w14, W_c=wc14),
                          "16": dict(W_star=w16, W_c=wc16)}
    d1s, d1c = w14 - w12, wc14 - wc12
    d2s, d2c = w16 - w14, wc16 - wc14
    cum = abs((w16 - w12) - (wc16 - wc12))
    same = (np.sign(d1s) == np.sign(d1c)) and (np.sign(d2s) == np.sign(d2c))
    if same and cum <= 0.5:
        verdict = (f"TRACKING-CONFIRMED at three sizes: W* drifts with "
                   f"W_c (steps {d1s:+.2f}/{d1c:+.2f} then "
                   f"{d2s:+.2f}/{d2c:+.2f}; cumulative gap {cum:.2f} "
                   f"<= 0.5) -- the sensitivity peak rides the moving "
                   f"edge of the crossover; the platform statement's "
                   f"registered follow-up closes AFFIRMATIVE.")
    elif abs(w16 - w12) <= 0.25 and abs(wc16 - wc12) >= 0.75:
        verdict = "PINNED: W* stays fixed while the crossover moves."
    else:
        verdict = (f"UNRESOLVED: steps W* {d1s:+.2f},{d2s:+.2f} vs W_c "
                   f"{d1c:+.2f},{d2c:+.2f} (cum gap {cum:.2f}).")
    md = ["# Q2c-full -- L = 16 third point (RESULTS)\n",
          f"dim {dim}; {N_REAL} realizations/W; evr-certified.\n",
          "| L | W* | W_c |", "|---|---|---|"]
    for Ls in ("12", "14", "16"):
        s = results["summary"][Ls]
        md.append(f"| {Ls} | {s['W_star']:.2f} | {s['W_c']:.2f} |")
    md.append(f"\n**Reading: {verdict}**")
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    results["verdict"] = verdict
    with open(os.path.join(HERE, "RESULTS_L16.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_l16.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md[-4:]))


if __name__ == "__main__":
    main()
