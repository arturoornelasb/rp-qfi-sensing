#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q4 -- operator-structure scope law (RUNNER). See README.md (frozen)."""
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(gammas=[1.4, 1.8], sizes=[512, 1024, 2048, 4096],
           n_real={512: 48, 1024: 24, 2048: 12, 4096: 6},
           window=0.25, fit_sizes=3, band=8, rank=4, seed=31, tol=0.15)
STRUCTS = ["diag", "banded", "lowrank", "dense"]


def grp(N, gamma, rng):
    W = rng.standard_normal((N, N))
    W = (W + W.T) / np.sqrt(2.0)
    H = N ** (-gamma / 2.0) * W
    H[np.diag_indices(N)] += rng.standard_normal(N)
    return H


def make_ops(N, rng, band, rank):
    ops = {}
    ops["diag"] = ("d", rng.standard_normal(N))
    B = rng.standard_normal((N, N))
    mask = np.abs(np.subtract.outer(np.arange(N), np.arange(N))) <= band
    Vb = np.where(mask, B, 0.0)
    ops["banded"] = ("m", (Vb + Vb.T) / np.sqrt(2 * band))
    vk = rng.standard_normal((N, rank))
    ops["lowrank"] = ("m", (vk @ vk.T) / np.sqrt(rank * N))
    G = rng.standard_normal((N, N))
    ops["dense"] = ("m", (G + G.T) / (np.sqrt(2.0) * np.sqrt(N)))
    return ops


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    cells = {}
    for gamma in cfg["gammas"]:
        for N in cfg["sizes"]:
            t0 = time.time()
            half = cfg["window"] / 2
            w0, w1 = int((0.5 - half) * N), int((0.5 + half) * N)
            cols = np.arange(w0, w1)
            acc = {s: [] for s in STRUCTS}
            ipr = []
            for _ in range(cfg["n_real"][N]):
                H = grp(N, gamma, rng)
                E, U = np.linalg.eigh(H)
                ipr.extend(np.log(np.sum(U[:, cols] ** 4, axis=0)).tolist())
                dE = E[:, None] - E[None, cols]
                np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
                inv2 = 1.0 / dE ** 2
                ops = make_ops(N, rng, cfg["band"], cfg["rank"])
                for s, (kind, V) in ops.items():
                    Mw = (U.T @ (V[:, None] * U[:, cols]) if kind == "d"
                          else U.T @ (V @ U[:, cols]))
                    acc[s].extend(np.log(np.sum(Mw ** 2 * inv2,
                                                axis=0)).tolist())
            cells[(gamma, N)] = dict(
                mln={s: float(np.mean(acc[s])) for s in STRUCTS},
                mln_ipr=float(np.mean(ipr)))
            print(f"[g={gamma} N={N}] " + " ".join(
                f"{s}:{cells[(gamma, N)]['mln'][s]:+.2f}" for s in STRUCTS)
                + f" ({time.time()-t0:.0f} s)", flush=True)

    md = ["# Q4 -- operator-structure scope law (RESULTS)\n",
          "| gamma | structure | alpha | enh = a - a_dense | 1 - D2 | "
          "class |", "|---|---|---|---|---|---|"]
    results = {"cells": {f"{g}|{N}": c for (g, N), c in cells.items()}}
    table = {}
    for gamma in cfg["gammas"]:
        Ns = cfg["sizes"][-cfg["fit_sizes"]:]
        x = np.log(Ns)
        d2 = float(-np.polyfit(x, [cells[(gamma, N)]["mln_ipr"]
                                   for N in Ns], 1)[0])
        alphas = {s: float(np.polyfit(
            x, [cells[(gamma, N)]["mln"][s] for N in Ns], 1)[0])
            for s in STRUCTS}
        for s in STRUCTS:
            enh = alphas[s] - alphas["dense"]
            if abs(enh - (1 - d2)) <= cfg["tol"]:
                cls = "INHERITS"
            elif abs(enh) < cfg["tol"]:
                cls = "BLIND"
            else:
                cls = f"PARTIAL ({enh/(1-d2):.0%})"
            table.setdefault(s, []).append(cls)
            md.append(f"| {gamma} | {s} | {alphas[s]:+.3f} | {enh:+.3f} | "
                      f"{1-d2:.3f} | {cls} |")
        results[f"fits_g{gamma}"] = dict(alphas=alphas, d2=d2)
    md.append("\n**Classification (consistent across gammas):**")
    for s in STRUCTS:
        labs = table[s]
        lab = labs[0] if all(
            l.split()[0] == labs[0].split()[0] for l in labs) \
            else f"MIXED {labs}"
        md.append(f"- {s}: **{lab}**")
        results.setdefault("classification", {})[s] = lab
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-8:]))


if __name__ == "__main__":
    main()
