#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q3c -- stratified gasket null (see README.md, frozen 2026-07-12).
Reuses the Q3b certified instrument; tags every sample by class
(reg/mult) and by energy window; per-stratum enhancement exponents."""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), "..", "q03b_certified_multiplets"))
from run_q03b import CFG as Q3BCFG, gasket_graph, certified_eigh, find_clusters

HERE = os.path.dirname(os.path.abspath(__file__))
E_CENTERS = [0.30, 0.43, 0.57, 0.70]
E_WIDTH = 0.10
CAP = 256
STRATA = ["reg", "mult"] + [f"E{c:.2f}" for c in E_CENTERS]


def collect(g, cfg, rng, results):
    nv, edges = gasket_graph(g)
    rows = [e[0] for e in edges] + [e[1] for e in edges]
    colz = [e[1] for e in edges] + [e[0] for e in edges]
    H = sp.coo_matrix((np.ones(len(rows)), (rows, colz)),
                      shape=(nv, nv)).toarray()
    Eref, U, cert = certified_eigh(H, cfg)
    rec = dict(N=nv, cert=cert)
    if not cert["ok"]:
        rec["status"] = "INSTRUMENT-FAIL"
        return rec
    N = nv
    tol_deg = cfg["clus_rel_tol"] * cert["width"]
    clusters = find_clusters(Eref, tol_deg)

    # stratum -> list of clusters (s, e)
    strata_cl = {s: [] for s in STRATA}
    half = cfg["window"] / 2.0
    w0, w1 = int((0.5 - half) * N), int((0.5 + half) * N)
    for s, e in clusters:
        mid = (s + e - 1) / 2
        if w0 <= mid < w1:
            strata_cl["reg" if e - s == 1 else "mult"].append((s, e))
        for c in E_CENTERS:
            a, b = int((c - E_WIDTH / 2) * N), int((c + E_WIDTH / 2) * N)
            if a <= mid < b:
                strata_cl[f"E{c:.2f}"].append((s, e))
    for k in strata_cl:
        picked, tot = [], 0
        for i in rng.permutation(len(strata_cl[k])):
            picked.append(strata_cl[k][i])
            tot += picked[-1][1] - picked[-1][0]
            if tot >= CAP:
                break
        strata_cl[k] = picked

    acc = {k: {"d": [], "g": [], "ipr": []} for k in STRATA}
    for _ in range(cfg["n_vdraws"]):
        vd = rng.standard_normal(N)
        Vg = rng.standard_normal((N, N))
        Vg = (Vg + Vg.T) / (np.sqrt(2.0) * np.sqrt(N))
        for k, picked in strata_cl.items():
            for s, e in picked:
                Uc = U[:, s:e]
                for tag, Mc in (("d", U.T @ (vd[:, None] * Uc)),
                                ("g", U.T @ (Vg @ Uc))):
                    Vcc = Mc[s:e, :]
                    _, R = np.linalg.eigh(0.5 * (Vcc + Vcc.T))
                    B = Mc @ R
                    dE = Eref[:, None] - Eref[s:e][None, :]
                    np.putmask(dE, np.abs(dE) < tol_deg, np.inf)
                    acc[k][tag].extend(
                        np.log(np.sum(B ** 2 / dE ** 2, axis=0)).tolist())
                    if tag == "d":
                        acc[k]["ipr"].extend(np.log(
                            np.sum((Uc @ R) ** 4, axis=0)).tolist())
        del Vg
    rec["status"] = "OK"
    rec["strata"] = {}
    for k in STRATA:
        if acc[k]["d"]:
            rec["strata"][k] = dict(
                n=len(acc[k]["d"]),
                mln_d=float(np.mean(acc[k]["d"])),
                mln_g=float(np.mean(acc[k]["g"])),
                mln_ipr=float(np.mean(acc[k]["ipr"])))
    return rec


def main():
    t00 = time.time()
    cfg = Q3BCFG
    rng = np.random.default_rng(21)
    results = {"sizes": {}}
    for g in cfg["generations"]:
        t0 = time.time()
        rec = collect(g, cfg, rng, results)
        results["sizes"][str(g)] = rec
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1, default=float)
        print(f"[g={g}] {rec.get('status')} ({time.time()-t0:.0f} s, "
              f"total {time.time()-t00:.0f} s)", flush=True)

    ok_g = [g for g in cfg["generations"]
            if results["sizes"][str(g)].get("status") == "OK"]
    fit_g = ok_g[-cfg["fit_sizes"]:]
    lnN = np.log([results["sizes"][str(g)]["N"] for g in fit_g])
    md = ["# Q3c -- stratified gasket null (RESULTS)\n",
          f"Sizes fitted: g = {fit_g}; strata: reg/mult + 4 energy "
          f"windows; cap {CAP} states/stratum.\n",
          "| stratum | enh (alpha_d - alpha_g) | D2_meas | dev vs 1-D2 |",
          "|---|---|---|---|"]
    results["fits"] = {}
    hidden = []
    for k in STRATA:
        try:
            ld = [results["sizes"][str(g)]["strata"][k]["mln_d"]
                  for g in fit_g]
            lg = [results["sizes"][str(g)]["strata"][k]["mln_g"]
                  for g in fit_g]
            li = [results["sizes"][str(g)]["strata"][k]["mln_ipr"]
                  for g in fit_g]
        except KeyError:
            md.append(f"| {k} | (insufficient) | - | - |")
            continue
        a_d = float(np.polyfit(lnN, ld, 1)[0])
        a_g = float(np.polyfit(lnN, lg, 1)[0])
        d2 = float(-np.polyfit(lnN, li, 1)[0])
        enh = a_d - a_g
        dev = abs(enh - (1.0 - d2))
        results["fits"][k] = dict(enh=enh, D2=d2, dev=dev)
        md.append(f"| {k} | {enh:+.3f} | {d2:.3f} | {dev:.3f} |")
        if abs(enh) >= cfg["r1_margin"]:
            hidden.append((k, enh))
    if not hidden:
        verdict = ("STRATIFIED-NULL-CONFIRMED: every stratum (class and "
                   "energy) has |enh| < 0.15 -- the Q3b null is not a "
                   "composition artifact; geometric fractality confers "
                   "no enhancement in any population. The paper's "
                   "criterion (randomness required) is strengthened.")
    else:
        verdict = ("HIDDEN-STRATUM: " + ", ".join(
            f"{k} enh {e:+.2f}" for k, e in hidden)
            + " -- see table (negative = Anderson-like suppressed "
              "branch).")
    if "reg" in results["fits"]:
        verdict += (f" R3: reg-stratum dev vs 1-D2 = "
                    f"{results['fits']['reg']['dev']:.3f} "
                    f"(Q3b pooled: 0.816).")
    md.append(f"\n**Reading: {verdict}**")
    md.append(f"\nWall: {time.time()-t00:.0f} s.")
    results["verdict"] = verdict
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1, default=float)
    print("\n".join(md[-5:]))


if __name__ == "__main__":
    main()
