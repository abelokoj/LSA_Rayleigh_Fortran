#!/usr/bin/env python3
"""Re-analyse the Gram matrix that the solver's LSA_activate==2 step prints, and compare two ways of choosing the most dangerous mode.

The solver (main.cpp, LSA_activate==2) prints one line per inner product of the Krylov members,
    i,j,dot_product <i> <j> <value>
and then fort_process_LSA (PROB_CPP_PARMS.F90) prints the eigenvalues of the companion matrix S,
    eigenvalue <i> = (<re>,<im>)
and its choice (max_index = the eigenvalue of largest real part). This script rebuilds G from the log and solves the same problem two ways:

  * as fort_process_LSA does: the pseudo-inverse of B^T W B keeping every positive singular value, and the eigenvalue of largest real part;
  * as svd_and_diagonalize_sussman.f90 does: the pseudo-inverse with a relative cutoff (rcond), and the growing, near-stationary mode of largest amplitude in v_1 (amplitude |b_j| * sqrt(x_j^H G x_j), X b = e_1).

It prints both, with s = ln(mu)/tau, against Chandrasekhar's exact viscous rate and the two-fluid viscous rate when --kr0 is given.

Member normalisation. The companion relation A V = V S holds for the unscaled sequence v_{k+1} = A v_k. If the solver rescales each member before it is re-added (in the version of 2026-10-15 the level set of each member is rescaled to a maximum of 2 dx in the band and the velocity to LSA_velocity_scale), the inner products are those of rescaled vectors and the eigenvalues of S are not those of A. With --scales FILE (one factor per member, true v_k = factor_k * stored v_k) the Gram matrix is corrected first; a single factor per member is only correct if every field was scaled by the same factor.

Usage:
    python3 krylov_test/lsa_log_spectrum.py lsa2.out --tau 0.25 --kr0 0.7 [--scales scales.txt] [--rcond 1e-12]
"""
import argparse
import re
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "plotting"))

DOT = re.compile(r"i,j,dot_product\s+(\d+)\s+(\d+)\s+([-+0-9.eEdD]+)")
EIG = re.compile(r"eigenvalue\s+(\d+)\s*=\s*\(\s*([-+0-9.eEdD]+)\s*,\s*([-+0-9.eEdD]+)\s*\)")


def parse_log(path):
    entries, eigs = {}, {}
    with open(path, errors="replace") as fh:
        for line in fh:
            m = DOT.search(line)
            if m:
                entries[(int(m.group(1)), int(m.group(2)))] = float(m.group(3).replace("D", "E").replace("d", "e"))
                continue
            m = EIG.search(line)
            if m:
                eigs[int(m.group(1))] = complex(float(m.group(2).replace("D", "E")), float(m.group(3).replace("D", "E")))
    if not entries:
        raise SystemExit("no 'i,j,dot_product' lines found in " + path)
    n = max(max(i, j) for i, j in entries) + 1
    G = np.full((n, n), np.nan)
    for (i, j), v in entries.items():
        G[i, j] = v
    if np.isnan(G).any():
        raise SystemExit("the Gram matrix in the log is incomplete")
    return G, [eigs[k] for k in sorted(eigs)]


def companion_spectrum(G, rcond, select):
    n = G.shape[0] - 1
    Gb, b = G[:n, :n], G[:n, n]
    w, V = np.linalg.eigh(0.5 * (Gb + Gb.T))      # SVD of a symmetric matrix
    sv = np.abs(w)
    keep = sv > (rcond * sv.max() if rcond > 0 else 0.0)
    a = V[:, keep] @ ((V[:, keep].T @ b) / w[keep])
    S = np.zeros((n, n))
    S[np.arange(1, n), np.arange(n - 1)] = 1.0
    S[:, -1] = a
    mu, X = np.linalg.eig(S)
    order = np.argsort(-np.abs(mu))
    mu, X = mu[order], X[:, order]
    coef = np.linalg.solve(X, np.eye(n)[:, 0])
    pn = np.sqrt(np.maximum(np.real(np.einsum("ij,ik,kj->j", X.conj(), Gb, X)), 0))
    amp = np.abs(coef) * pn
    resid = np.sqrt(max(G[n, n] - 2 * a @ b + a @ Gb @ a, 0) / G[n, n]) if G[n, n] > 0 else np.nan
    return mu, amp, int(keep.sum()), resid


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("log")
    ap.add_argument("--tau", type=float, required=True, help="time between members (fixed_dt * steps per member)")
    ap.add_argument("--kr0", type=float, default=None)
    ap.add_argument("--mu", type=float, default=0.02)
    ap.add_argument("--rho-g", type=float, default=1.225e-3, dest="rho_g")
    ap.add_argument("--mu-g", type=float, default=2.6e-4, dest="mu_g")
    ap.add_argument("--scales", default=None, help="file with one factor per member: true v_k = factor_k * stored v_k")
    ap.add_argument("--rcond", type=float, default=1e-12)
    a = ap.parse_args()

    G, printed = parse_log(a.log)
    N = G.shape[0]
    if a.scales:
        f = np.loadtxt(a.scales, ndmin=1)
        if f.size != N:
            raise SystemExit(f"--scales has {f.size} factors, the log has {N} members")
        G = G * np.outer(f, f)
    print(f"{N} Krylov members, tau = {a.tau}")

    ref = {}
    if a.kr0:
        import dispersion
        from plot_mode_shape_archive import growth_rate_two_fluid
        ref["Chandrasekhar"] = float(dispersion.growth_rate_viscous_exact(a.kr0, mu=a.mu))
        ref["two-fluid viscous"] = growth_rate_two_fluid(a.kr0, a.mu, a.rho_g, a.mu_g / a.mu)

    def report(title, mu, amp, kept, resid, pick):
        s = np.log(mu.astype(complex)) / a.tau
        print(f"\n{title}: rank {kept} of {N - 1}, residual {resid:.3e}")
        print("   j        |mu|      Re s        Im s     amplitude")
        for j in range(min(len(mu), 12)):
            mark = "  <-- selected" if j == pick else ""
            print(f"{j + 1:4d} {abs(mu[j]):11.6f} {s[j].real:9.5f} {s[j].imag:11.5f} {amp[j]:12.4e}{mark}")
        if pick is not None:
            for k, v in ref.items():
                print(f"   selected s = {s[pick].real:.6f}  vs {k} {v:.6f}: {100 * (s[pick].real / v - 1):+.2f} %")

    mu1, amp1, k1, r1 = companion_spectrum(G, 0.0, "maxre")
    pick1 = int(np.argmax(mu1.real))
    report("as fort_process_LSA (no cutoff, largest Re mu)", mu1, amp1, k1, r1, pick1)
    mu2, amp2, k2, r2 = companion_spectrum(G, a.rcond, "dominant")
    s2 = np.log(mu2.astype(complex)) / a.tau
    ok = (s2.real > 0) & (np.abs(s2.imag) <= 0.1 * s2.real)
    pick2 = int(np.argmax(np.where(ok, amp2, -1))) if ok.any() else None
    report(f"as svd_and_diagonalize_sussman (rcond {a.rcond:g}, dominant amplitude)", mu2, amp2, k2, r2, pick2)
    if printed:
        print(f"\nthe solver printed {len(printed)} eigenvalues; largest Re mu = {max(printed, key=lambda z: z.real):.6f}")


if __name__ == "__main__":
    main()
