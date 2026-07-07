#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q2 / DR-A -- platform coexistence test (RUNNER).
See README.md (preregistered). t-V chain + disorder; tilt = sensed
parameter at h = 0; chi_f per mid-window eigenstate; r-tilde alongside."""
import json
import os
import time
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    J=1.0, V=1.0,
    W_grid=[0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0],
    sizes={12: 200, 14: 64},
    window=0.25,
    tail_frac=0.05,
    cert_tol=1e-12,
    seed=21,
    r1_lo=0.42, r1_hi=0.50, r1_erg=0.52, r1_loc=0.40,
    r2_lo=1.6, r2_hi=2.4,
)


def basis_half_filling(L):
    states = [sum(1 << i for i in occ)
              for occ in combinations(range(L), L // 2)]
    states.sort()
    index = {s: k for k, s in enumerate(states)}
    return states, index


def build_h(L, states, index, w, J, V):
    dim = len(states)
    H = np.zeros((dim, dim))
    for k, s in enumerate(states):
        diag = 0.0
        for i in range(L):
            if s >> i & 1:
                diag += w[i]
                if i < L - 1 and (s >> (i + 1) & 1):
                    diag += V
        H[k, k] = diag
        for i in range(L - 1):
            if (s >> i & 1) and not (s >> (i + 1) & 1):
                s2 = s ^ (1 << i) ^ (1 << (i + 1))
                k2 = index[s2]
                H[k2, k] += -J
                H[k, k2] += -J
    return H


def dipole_diag(L, states):
    c = (L - 1) / 2.0
    d = np.zeros(len(states))
    for k, s in enumerate(states):
        d[k] = sum((i - c) for i in range(L) if s >> i & 1)
    return d


def rt_from_levels(E):
    s = np.diff(E)
    s = s[s > 0]
    if len(s) < 2:
        return np.empty(0)
    return np.minimum(s[:-1], s[1:]) / np.maximum(s[:-1], s[1:])


def hill_index(x, frac):
    x = np.sort(np.asarray(x))
    k = max(20, int(frac * len(x)))
    tail = x[-k:]
    return 1.0 + 1.0 / float(np.mean(np.log(tail / tail[0]))), int(k)


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": {k: (v if not isinstance(v, dict) else
                              {str(a): b for a, b in v.items()})
                          for k, v in cfg.items()}, "cells": {}}

    for L, n_real in cfg["sizes"].items():
        states, index = basis_half_filling(L)
        dim = len(states)
        vop = dipole_diag(L, states)
        half = cfg["window"] / 2.0
        w0, w1 = int((0.5 - half) * dim), int((0.5 + half) * dim)
        cols = np.arange(w0, w1)
        for W in cfg["W_grid"]:
            t0 = time.time()
            rts, ln_chi, chis = [], [], []
            cert = None
            for rr in range(n_real):
                w = rng.uniform(-W, W, L)
                H = build_h(L, states, index, w, cfg["J"], cfg["V"])
                E, U = np.linalg.eigh(H)
                if rr == 0:
                    Ev = np.linalg.eigvalsh(H)
                    width = float(Ev[-1] - Ev[0])
                    agree = float(np.max(np.abs(E - Ev)))
                    cert = dict(agree=agree, width=width,
                                ok=bool(agree <= cfg["cert_tol"] * width))
                rts.extend(rt_from_levels(E[w0:w1]).tolist())
                dE = E[:, None] - E[None, cols]
                np.putmask(dE, np.abs(dE) < 1e-14, np.inf)
                Mw = U.T @ (vop[:, None] * U[:, cols])
                chi = np.sum(Mw ** 2 / dE ** 2, axis=0)
                ln_chi.extend(np.log(chi).tolist())
                chis.extend(chi.tolist())
            rts = np.array(rts)
            aH, ktail = hill_index(chis, cfg["tail_frac"])
            cell = dict(dim=dim, n_real=n_real, cert=cert,
                        r_mean=float(rts.mean()),
                        r_sem=float(rts.std() / np.sqrt(len(rts))),
                        n_ratios=int(len(rts)),
                        mln_chi=float(np.mean(ln_chi)),
                        mean_chi=float(np.mean(chis)),
                        hill=float(aH), hill_k=ktail)
            results["cells"][f"{L}|{W}"] = cell
            print(f"[L={L} W={W:>3}] <r> = {cell['r_mean']:.4f}"
                  f"({cell['r_sem']:.4f}), ln chi_typ = "
                  f"{cell['mln_chi']:+.3f}, Hill = {aH:.2f} "
                  f"(cert {cert['agree']:.1e}, {time.time()-t0:.1f} s)")
            with open(os.path.join(HERE, "results_raw.json"), "w") as f:
                json.dump(results, f, indent=1)

    # ---------------- preregistered verdicts (L = 14) ----------------
    L14 = {float(k.split("|")[1]): results["cells"][k]
           for k in results["cells"] if k.startswith("14|")}
    Ws = sorted(L14)
    Wstar = max(Ws, key=lambda W: L14[W]["mln_chi"])
    rstar = L14[Wstar]["r_mean"]
    if cfg["r1_lo"] <= rstar <= cfg["r1_hi"]:
        r1 = f"SUPPORTS (W* = {Wstar:g}, <r> = {rstar:.4f} mid-crossover)"
    elif rstar >= cfg["r1_erg"]:
        r1 = f"CHALLENGES-ERGODIC (W* = {Wstar:g}, <r> = {rstar:.4f})"
    elif rstar <= cfg["r1_loc"]:
        r1 = f"CHALLENGES-LOCALIZED (W* = {Wstar:g}, <r> = {rstar:.4f})"
    else:
        r1 = (f"BOUNDARY (W* = {Wstar:g}, <r> = {rstar:.4f} between the "
              f"gated bands)")
    aH_star = L14[Wstar]["hill"]
    r2 = (f"SUPPORTS (Hill = {aH_star:.2f} at W*)"
          if cfg["r2_lo"] <= aH_star <= cfg["r2_hi"]
          else f"outside band (Hill = {aH_star:.2f} at W*)")

    md = ["# Q2 / DR-A -- platform coexistence test (RESULTS)\n",
          f"t-V chain, J = V = 1, half filling; L = 12 (924, 200 real) and "
          f"L = 14 (3432, 64 real); window central 25%; chi_f w.r.t. the "
          f"tilt at h = 0.\n",
          "| L | W | <r> | ln chi_typ | mean chi | Hill (top 5%) | cert |",
          "|---|---|---|---|---|---|---|"]
    for key in sorted(results["cells"],
                      key=lambda k: (int(k.split("|")[0]),
                                     float(k.split("|")[1]))):
        c = results["cells"][key]
        Lk, Wk = key.split("|")
        md.append(f"| {Lk} | {Wk} | {c['r_mean']:.4f}({c['r_sem']:.4f}) | "
                  f"{c['mln_chi']:+.3f} | {c['mean_chi']:.3g} | "
                  f"{c['hill']:.2f} | {c['cert']['agree']:.0e} |")
    md.append(f"\n- **R1 (coexistence)**: {r1}")
    md.append(f"- **R2 (GOE x^-2 tail)**: {r2}")
    md.append(f"- Ergodic control (W = 0.5, L = 14): <r> = "
              f"{L14[0.5]['r_mean']:.4f}, Hill = {L14[0.5]['hill']:.2f}")
    verdict = ("SUPPORTS the platform minimal model"
               if r1.startswith("SUPPORTS") and r2.startswith("SUPPORTS")
               else ("PARTIAL -- see R1/R2"
                     if not (r1.startswith("CHALLENGES")) else
                     "CHALLENGES the platform minimal model"))
    md.append(f"\n**Reading: {verdict}.**")
    results["verdict"] = dict(R1=r1, R2=r2, overall=verdict,
                              W_star=float(Wstar))
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
