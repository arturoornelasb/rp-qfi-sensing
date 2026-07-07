#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q3 / DR-C -- QFI on the Sierpinski gasket (RUNNER).
See README.md (preregistered). Fractality from geometry, zero disorder;
statistics from a V-ensemble on the fixed lattice."""
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    generations=[5, 6, 7, 8],
    fit_sizes=3,
    n_vdraws=8,
    window=0.25,
    max_chi_states=512,
    min_usable=64,
    degen_rel_tol=1e-10,
    seed=11,
    r1_margin=0.15,
    r2_tol=0.15,
)


def gasket_graph(g):
    """Vertices and edges of the level-g Sierpinski gasket."""
    h = np.sqrt(3.0) / 2.0
    tris = [((0.0, 0.0), (1.0, 0.0), (0.5, h))]
    for _ in range(g):
        nxt = []
        for a, b, c in tris:
            ab = ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
            bc = ((b[0] + c[0]) / 2, (b[1] + c[1]) / 2)
            ca = ((c[0] + a[0]) / 2, (c[1] + a[1]) / 2)
            nxt += [(a, ab, ca), (ab, b, bc), (ca, bc, c)]
        tris = nxt
    idx, edges = {}, set()
    def vid(p):
        key = (round(p[0], 9), round(p[1], 9))
        if key not in idx:
            idx[key] = len(idx)
        return idx[key]
    for a, b, c in tris:
        ia, ib, ic = vid(a), vid(b), vid(c)
        edges |= {tuple(sorted(e)) for e in [(ia, ib), (ib, ic), (ic, ia)]}
    return len(idx), sorted(edges)


def main():
    t00 = time.time()
    cfg = CFG
    rng = np.random.default_rng(cfg["seed"])
    results = {"config": dict(cfg), "sizes": {}}

    for g in cfg["generations"]:
        t0 = time.time()
        N, edges = gasket_graph(g)
        H = np.zeros((N, N))
        for i, j in edges:
            H[i, j] = H[j, i] = -1.0
        E, U = np.linalg.eigh(H)
        width = float(E[-1] - E[0])
        half = cfg["window"] / 2.0
        w0, w1 = int((0.5 - half) * N), int((0.5 + half) * N)
        gaps = np.minimum(np.abs(np.diff(E, prepend=E[0] - 1.0)),
                          np.abs(np.diff(E, append=E[-1] + 1.0)))
        usable = np.arange(w0, w1)[gaps[w0:w1] > cfg["degen_rel_tol"] * width]
        degen_frac = 1.0 - len(usable) / (w1 - w0)
        if len(usable) > cfg["max_chi_states"]:
            usable = np.sort(rng.choice(usable, cfg["max_chi_states"],
                                        replace=False))
        ipr = np.sum(U[:, usable] ** 4, axis=0)
        ln_chi_d, ln_chi_g = [], []
        for _ in range(cfg["n_vdraws"]):
            vd = rng.standard_normal(N)
            Vg = rng.standard_normal((N, N))
            Vg = (Vg + Vg.T) / (np.sqrt(2.0) * np.sqrt(N))
            dE = E[:, None] - E[None, usable]
            np.putmask(dE, np.abs(dE) < cfg["degen_rel_tol"] * width, np.inf)
            Md = U.T @ (vd[:, None] * U[:, usable])
            Mg = U.T @ (Vg @ U[:, usable])
            ln_chi_d.extend(np.log(np.sum(Md ** 2 / dE ** 2, axis=0)).tolist())
            ln_chi_g.extend(np.log(np.sum(Mg ** 2 / dE ** 2, axis=0)).tolist())
            del Vg, Md, Mg, dE
        results["sizes"][str(g)] = dict(
            N=N, n_edges=len(edges), width=width,
            n_window=int(w1 - w0), n_usable=int(len(usable)),
            degen_frac=float(degen_frac),
            mln_ipr=float(np.mean(np.log(ipr))),
            mln_chi_d=float(np.mean(ln_chi_d)),
            mln_chi_g=float(np.mean(ln_chi_g)),
            limited=bool(len(usable) < cfg["min_usable"]))
        print(f"[g={g} N={N:>5}] usable {len(usable)}/{w1-w0} "
              f"(degen frac {degen_frac:.2f}), ln chi_d "
              f"{results['sizes'][str(g)]['mln_chi_d']:+.3f} "
              f"({time.time()-t0:.1f} s)")
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1)

    # ---------------- fits + preregistered verdicts ----------------
    gs = cfg["generations"][-cfg["fit_sizes"]:]
    x = np.array([np.log(results["sizes"][str(g)]["N"]) for g in gs])
    yd = np.array([results["sizes"][str(g)]["mln_chi_d"] for g in gs])
    yg = np.array([results["sizes"][str(g)]["mln_chi_g"] for g in gs])
    yi = np.array([results["sizes"][str(g)]["mln_ipr"] for g in gs])
    alpha_d = float(np.polyfit(x, yd, 1)[0])
    alpha_g = float(np.polyfit(x, yg, 1)[0])
    d2 = float(-np.polyfit(x, yi, 1)[0])
    enh = alpha_d - alpha_g
    limited = any(results["sizes"][str(g)]["limited"] for g in gs)

    if limited:
        r1 = "DEGENERACY-LIMITED (usable states below floor at a fit size)"
    elif enh >= cfg["r1_margin"]:
        r1 = "ENHANCED"
    elif enh <= -cfg["r1_margin"]:
        r1 = "SUPPRESSED (Anderson-like sign)"
    else:
        r1 = "NO EFFECT"
    dev = abs(enh - (1.0 - d2))
    r2 = ("RP FORMULA EXTENDS" if dev <= cfg["r2_tol"] and not limited
          else "DOES NOT EXTEND (report the measured pair)")

    md = ["# Q3 / DR-C -- Sierpinski gasket QFI (RESULTS)\n",
          f"Geometry-only fractality, V-ensemble ({cfg['n_vdraws']} draws), "
          f"fits over g = {gs} (largest {cfg['fit_sizes']} sizes).\n",
          "| g | N | usable/window | degen frac | ln chi_d | ln chi_g | "
          "ln IPR |", "|---|---|---|---|---|---|---|"]
    for g in cfg["generations"]:
        s = results["sizes"][str(g)]
        md.append(f"| {g} | {s['N']} | {s['n_usable']}/{s['n_window']} | "
                  f"{s['degen_frac']:.2f} | {s['mln_chi_d']:+.3f} | "
                  f"{s['mln_chi_g']:+.3f} | {s['mln_ipr']:+.3f} |")
    md.append(f"\n- alpha_diag = {alpha_d:+.3f}; alpha_goe = {alpha_g:+.3f}; "
              f"enh = {enh:+.3f}; D2_meas = {d2:.3f}; 1 - D2_meas = "
              f"{1-d2:.3f}; dev = {dev:.3f}")
    md.append(f"- **R1**: {r1} (margin {cfg['r1_margin']})")
    md.append(f"- **R2**: {r2} (tol {cfg['r2_tol']})")
    md.append(f"\n**Reading: {r1}; {r2}.**")
    results["fits"] = dict(alpha_diag=alpha_d, alpha_goe=alpha_g, enh=enh,
                           d2_meas=d2, dev=dev)
    results["verdict"] = dict(R1=r1, R2=r2)
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
