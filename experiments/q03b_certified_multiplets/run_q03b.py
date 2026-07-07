#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Q3b -- Sierpinski gasket QFI, certified instrument + degenerate-multiplet
subspace treatment (RUNNER). See README.md (preregistered)."""
import json
import os
import time

import numpy as np
from scipy.linalg import eigh as seigh

HERE = os.path.dirname(os.path.abspath(__file__))

CFG = dict(
    generations=[5, 6, 7, 8],
    fit_sizes=3,
    n_vdraws=8,
    window=0.25,
    max_chi_states=512,
    clus_rel_tol=1e-10,
    cert_val_tol=1e-12,
    cert_res_tol=1e-10,
    seed=12,
    r1_margin=0.15,
    r2_tol=0.15,
)


def gasket_graph(g):
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


def certified_eigh(H, cfg):
    Eref = np.linalg.eigvalsh(H)
    E, U = seigh(H, driver="evr")
    width = float(Eref[-1] - Eref[0])
    agree = float(np.max(np.abs(E - Eref)))
    cols = np.linspace(0, len(E) - 1, 64).astype(int)
    res = float(np.max(np.abs(H @ U[:, cols] - U[:, cols] * E[cols])))
    ok = (agree <= cfg["cert_val_tol"] * width
          and res <= cfg["cert_res_tol"] * width)
    return Eref, U, dict(width=width, agree=agree, resid=res, ok=bool(ok))


def find_clusters(E, tol_abs):
    br = np.where(np.diff(E) > tol_abs)[0]
    starts = np.r_[0, br + 1]
    ends = np.r_[br + 1, len(E)]
    return list(zip(starts.tolist(), ends.tolist()))


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
        Eref, U, cert = certified_eigh(H, cfg)
        rec = dict(N=N, cert=cert)
        if not cert["ok"]:
            rec["status"] = "INSTRUMENT-FAIL"
            results["sizes"][str(g)] = rec
            print(f"[g={g} N={N}] INSTRUMENT-FAIL "
                  f"(agree {cert['agree']:.1e}, resid {cert['resid']:.1e})")
            continue
        tol_deg = cfg["clus_rel_tol"] * cert["width"]
        half = cfg["window"] / 2.0
        w0, w1 = int((0.5 - half) * N), int((0.5 + half) * N)
        wcl = [c for c in find_clusters(Eref, tol_deg)
               if w0 <= (c[0] + c[1] - 1) / 2 < w1]
        picked, tot = [], 0
        for i in rng.permutation(len(wcl)):
            s, e = wcl[i]
            picked.append((int(s), int(e)))
            tot += e - s
            if tot >= cfg["max_chi_states"]:
                break
        ln_chi = {"d": [], "g": []}
        ln_ipr, ln_ipr_nd = [], []
        for _ in range(cfg["n_vdraws"]):
            vd = rng.standard_normal(N)
            Vg = rng.standard_normal((N, N))
            Vg = (Vg + Vg.T) / (np.sqrt(2.0) * np.sqrt(N))
            for s, e in picked:
                Uc = U[:, s:e]
                for tag, Mc in (("d", U.T @ (vd[:, None] * Uc)),
                                ("g", U.T @ (Vg @ Uc))):
                    Vcc = Mc[s:e, :]
                    _, R = np.linalg.eigh(0.5 * (Vcc + Vcc.T))
                    B = Mc @ R
                    dE = Eref[:, None] - Eref[s:e][None, :]
                    np.putmask(dE, np.abs(dE) < tol_deg, np.inf)
                    ln_chi[tag].extend(
                        np.log(np.sum(B ** 2 / dE ** 2, axis=0)).tolist())
                    if tag == "d":
                        ln_ipr.extend(np.log(
                            np.sum((Uc @ R) ** 4, axis=0)).tolist())
                        if e - s == 1:
                            ln_ipr_nd.append(float(np.log(
                                np.sum(Uc[:, 0] ** 4))))
            del Vg
        rec.update(status="OK", n_states=int(tot),
                   n_multiplets=len(picked),
                   max_multiplet=int(max(e - s for s, e in picked)),
                   mln_chi_d=float(np.mean(ln_chi["d"])),
                   mln_chi_g=float(np.mean(ln_chi["g"])),
                   mln_ipr=float(np.mean(ln_ipr)),
                   mln_ipr_nondeg=(float(np.mean(ln_ipr_nd))
                                   if ln_ipr_nd else None))
        results["sizes"][str(g)] = rec
        print(f"[g={g} N={N:>5}] {tot} states / {len(picked)} multiplets "
              f"(max {rec['max_multiplet']}), ln chi_d "
              f"{rec['mln_chi_d']:+.3f}, ln IPR {rec['mln_ipr']:+.3f} "
              f"(agree {cert['agree']:.1e}, {time.time()-t0:.1f} s)")
        with open(os.path.join(HERE, "results_raw.json"), "w") as f:
            json.dump(results, f, indent=1)

    # ---------------- fits + preregistered verdicts ----------------
    valid = [g for g in cfg["generations"]
             if results["sizes"][str(g)].get("status") == "OK"]
    gs = valid[-cfg["fit_sizes"]:]
    coverage_ok = len(gs) >= cfg["fit_sizes"]
    md = ["# Q3b -- certified multiplet instrument (RESULTS)\n",
          f"Valid sizes: {valid} (fit over {gs}); V-ensemble "
          f"{cfg['n_vdraws']} draws; full-window multiplet treatment.\n",
          "| g | N | states/multiplets (max) | ln chi_d | ln chi_g | "
          "ln IPR | cert agree |", "|---|---|---|---|---|---|---|"]
    for g in cfg["generations"]:
        s = results["sizes"][str(g)]
        if s.get("status") != "OK":
            md.append(f"| {g} | {s['N']} | INSTRUMENT-FAIL | | | | "
                      f"{s['cert']['agree']:.1e} |")
            continue
        md.append(f"| {g} | {s['N']} | {s['n_states']}/{s['n_multiplets']} "
                  f"({s['max_multiplet']}) | {s['mln_chi_d']:+.3f} | "
                  f"{s['mln_chi_g']:+.3f} | {s['mln_ipr']:+.3f} | "
                  f"{s['cert']['agree']:.1e} |")
    if coverage_ok:
        x = np.array([np.log(results["sizes"][str(g)]["N"]) for g in gs])
        yd = np.array([results["sizes"][str(g)]["mln_chi_d"] for g in gs])
        yg = np.array([results["sizes"][str(g)]["mln_chi_g"] for g in gs])
        yi = np.array([results["sizes"][str(g)]["mln_ipr"] for g in gs])
        alpha_d = float(np.polyfit(x, yd, 1)[0])
        alpha_g = float(np.polyfit(x, yg, 1)[0])
        d2 = float(-np.polyfit(x, yi, 1)[0])
        enh = alpha_d - alpha_g
        dev = abs(enh - (1.0 - d2))
        if enh >= cfg["r1_margin"]:
            r1 = "ENHANCED"
        elif enh <= -cfg["r1_margin"]:
            r1 = "SUPPRESSED (Anderson-like sign)"
        else:
            r1 = "NO EFFECT"
        r2 = ("RP FORMULA EXTENDS" if dev <= cfg["r2_tol"]
              else "DOES NOT EXTEND (report the measured pair)")
        md.append(f"\n- alpha_diag = {alpha_d:+.3f}; alpha_goe = "
                  f"{alpha_g:+.3f}; enh = {enh:+.3f}; D2_meas = {d2:.3f}; "
                  f"1 - D2_meas = {1-d2:.3f}; dev = {dev:.3f}")
        md.append(f"- **R1**: {r1} (margin {cfg['r1_margin']})")
        md.append(f"- **R2**: {r2} (tol {cfg['r2_tol']})")
        md.append(f"\n**Reading: {r1}; {r2}.**")
        results["fits"] = dict(alpha_diag=alpha_d, alpha_goe=alpha_g,
                               enh=enh, d2_meas=d2, dev=dev)
        results["verdict"] = dict(R1=r1, R2=r2)
    else:
        md.append(f"\n**Reading: COVERAGE-LIMITED -- only {len(gs)} valid "
                  f"sizes (< {cfg['fit_sizes']}); no fits.**")
        results["verdict"] = dict(R1="COVERAGE-LIMITED", R2="COVERAGE-LIMITED")
    results["wall_time_s"] = round(time.time() - t00, 1)
    md.append(f"\nWall time: {results['wall_time_s']} s.")
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(HERE, "results_raw.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\n".join(md))


if __name__ == "__main__":
    main()
