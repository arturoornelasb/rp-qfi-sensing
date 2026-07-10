#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q5 -- PBRM third family (RUNNER). See README.md (frozen)."""
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = dict(bands=[0.5, 1.5, 4.0], sizes=[512, 1024, 2048, 4096],
           n_real={512: 48, 1024: 24, 2048: 12, 4096: 6},
           window=0.25, fit_sizes=3, seed=41, tol=0.15)


def pbrm(N, b, rng):
    i = np.arange(N)
    sig = 1.0 / np.sqrt(1.0 + (np.abs(np.subtract.outer(i, i)) / b) ** 2)
    W = rng.standard_normal((N, N))
    H = W * sig
    return (H + H.T) / np.sqrt(2.0)


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    cells = {}
    for b in cfg["bands"]:
        for N in cfg["sizes"]:
            t0 = time.time()
            half = cfg["window"] / 2
            w0, w1 = int((0.5 - half) * N), int((0.5 + half) * N)
            cols = np.arange(w0, w1)
            ln_d, ln_g, ipr = [], [], []
            for _ in range(cfg["n_real"][N]):
                H = pbrm(N, b, rng)
                E, U = np.linalg.eigh(H)
                ipr.extend(np.log(np.sum(U[:, cols] ** 4, axis=0)).tolist())
                dE = E[:, None] - E[None, cols]
                np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
                inv2 = 1.0 / dE ** 2
                vd = rng.standard_normal(N)
                G = rng.standard_normal((N, N))
                Vg = (G + G.T) / (np.sqrt(2.0) * np.sqrt(N))
                Md = U.T @ (vd[:, None] * U[:, cols])
                Mg = U.T @ (Vg @ U[:, cols])
                ln_d.extend(np.log(np.sum(Md ** 2 * inv2, axis=0)).tolist())
                ln_g.extend(np.log(np.sum(Mg ** 2 * inv2, axis=0)).tolist())
            cells[(b, N)] = dict(d=float(np.mean(ln_d)),
                                 g=float(np.mean(ln_g)),
                                 ipr=float(np.mean(ipr)))
            print(f"[b={b} N={N}] d {cells[(b, N)]['d']:+.2f} "
                  f"({time.time()-t0:.0f} s)", flush=True)

    md = ["# Q5 -- PBRM third family (RESULTS)\n",
          "| b | alpha_diag | alpha_dense | enh | D2_meas | 1 - D2 | "
          "dev |", "|---|---|---|---|---|---|---|"]
    results = {"cells": {f"{b}|{N}": c for (b, N), c in cells.items()}}
    ok = 0
    for b in cfg["bands"]:
        Ns = cfg["sizes"][-cfg["fit_sizes"]:]
        x = np.log(Ns)
        ad = float(np.polyfit(x, [cells[(b, N)]["d"] for N in Ns], 1)[0])
        ag = float(np.polyfit(x, [cells[(b, N)]["g"] for N in Ns], 1)[0])
        d2 = float(-np.polyfit(x, [cells[(b, N)]["ipr"] for N in Ns],
                               1)[0])
        enh = ad - ag
        dev = abs(enh - (1 - d2))
        ok += dev <= cfg["tol"]
        results[f"fit_b{b}"] = dict(alpha_diag=ad, alpha_dense=ag,
                                    d2=d2, enh=enh, dev=dev)
        md.append(f"| {b} | {ad:+.3f} | {ag:+.3f} | {enh:+.3f} | "
                  f"{d2:.3f} | {1-d2:.3f} | {dev:.3f} |")
    verdict = (f"EXTENDS to PBRM ({ok}/3 within tol): the 1 - D2 law is "
               f"about multifractal support, not RP minibands specifically."
               if ok >= 2 else
               f"DOES NOT EXTEND ({ok}/3 within tol): the enhancement is "
               f"RP-miniband-specific -- with Q1b (RP yes) and Q3b (gasket "
               f"no), the criterion is the RP structure, sharpened.")
    md.append(f"\n**Reading: {verdict}**")
    results["verdict"] = verdict
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
