"""
confined_dispersion.py
======================
Linear growth rate of the capillary (Rayleigh-Plateau) instability of a viscous liquid column of radius R0 inside a coaxial no-slip wall at r = Rw, with a viscous ambient fluid filling R0 < r < Rw. Axisymmetric perturbations proportional to exp(i k z + s t); capillary units R0 = sigma = rho_1 = 1, so the liquid viscosity is the Ohnesorge number.

Each fluid carries a potential part (pressure p = -rho s Phi) and a vortical part with l^2 = k^2 + s/nu:

    liquid:  u_r = A1 k I1(kr) + C1 I1(l1 r)
             u_z = i [A1 k I0(kr) + C1 (l1/k) I0(l1 r)]
             p   = -rho1 s A1 I0(kr)
    ambient: u_r = -A2 k K1(kr) + B2 k I1(kr) + C2 K1(l2 r) + D2 I1(l2 r)
             u_z = i [A2 k K0(kr) + B2 k I0(kr) - C2 (l2/k) K0(l2 r) + D2 (l2/k) I0(l2 r)]
             p   = -rho2 s [A2 K0(kr) + B2 I0(kr)]

Conditions at r = 1: continuity of u_r and u_z; continuity of tangential stress mu (du_z/dr + du_r/dz); normal-stress balance (-p + 2 mu du_r/dr)_1 - (-p + 2 mu du_r/dr)_2 = sigma (1 - k^2) eta, with the kinematic condition eta = u_r(1)/s. At r = Rw: u_r = u_z = 0 (no slip). The growth rate is the real root s of det M(s) = 0 (6 x 6); without the wall (B2 = D2 = 0, wall rows dropped) the 4 x 4 system is the unconfined two-fluid relation of Tomotika (1935).

Double precision (numpy/scipy), with exponentially scaled Bessel functions so that every entry stays finite for any viscosity (see _rows); each column of M is also normalised by its largest entry before the determinant is formed, so that the sign test that brackets the root is well scaled.
"""
import numpy as np
from scipy.special import iv, kv, ive, kve
from scipy.optimize import brentq

def _rows(s, k, oh, rho2, mu2, Rw):
    """6 x 6 (or 4 x 4 without the wall) boundary-condition matrix, columns A1, C1, A2, B2, C2, D2.

    Each basis function is multiplied by a constant chosen so that it is of order one where it is largest, which keeps every entry finite even when l*r is in the thousands (very small viscosity): I-type functions of a*r are scaled by exp(-a*R_ref) and K-type functions by exp(+a*1), where R_ref is 1 inside the liquid and R_w (or 1 without a wall) in the ambient. Scaling a column by a constant leaves the roots of det M unchanged. Values use scipy's exponentially scaled ive/kve: I_n(z) exp(-c) = ive(n, z) exp(z - c), K_n(z) exp(c) = kve(n, z) exp(c - z).
    """
    s, k = float(s), float(k)
    mu1, rho1 = float(oh), 1.0
    mu2, rho2 = float(mu2), float(rho2)
    nu1, nu2 = mu1 / rho1, mu2 / rho2
    l1 = np.sqrt(k * k + s / nu1); l2 = np.sqrt(k * k + s / nu2)
    Rref = 1.0 if Rw is None else float(Rw)
    def I(n, a, r, c):  return ive(n, a * r) * np.exp(a * r - c)
    def K(n, a, r, c):  return kve(n, a * r) * np.exp(c - a * r)
    # each basis: returns (u_r, W = u_z/i, dW/dr, du_r/dr, -p) at r, already scaled
    def liq_A1(r):
        c = k * 1.0
        return (k * I(1, k, r, c), k * I(0, k, r, c), k * k * I(1, k, r, c),
                k * k * (I(0, k, r, c) - I(1, k, r, c) / (k * r)), rho1 * s * I(0, k, r, c))
    def liq_C1(r):
        c = l1 * 1.0
        return (I(1, l1, r, c), (l1 / k) * I(0, l1, r, c), (l1 / k) * l1 * I(1, l1, r, c),
                l1 * (I(0, l1, r, c) - I(1, l1, r, c) / (l1 * r)), 0.0)
    def gas_A2(r):
        c = k * 1.0
        return (-k * K(1, k, r, c), k * K(0, k, r, c), -k * k * K(1, k, r, c),
                k * k * (K(0, k, r, c) + K(1, k, r, c) / (k * r)), rho2 * s * K(0, k, r, c))
    def gas_B2(r):
        c = k * Rref
        return (k * I(1, k, r, c), k * I(0, k, r, c), k * k * I(1, k, r, c),
                k * k * (I(0, k, r, c) - I(1, k, r, c) / (k * r)), rho2 * s * I(0, k, r, c))
    def gas_C2(r):
        c = l2 * 1.0
        return (K(1, l2, r, c), -(l2 / k) * K(0, l2, r, c), (l2 / k) * l2 * K(1, l2, r, c),
                -l2 * (K(0, l2, r, c) + K(1, l2, r, c) / (l2 * r)), 0.0)
    def gas_D2(r):
        c = l2 * Rref
        return (I(1, l2, r, c), (l2 / k) * I(0, l2, r, c), (l2 / k) * l2 * I(1, l2, r, c),
                l2 * (I(0, l2, r, c) - I(1, l2, r, c) / (l2 * r)), 0.0)
    gas = [gas_A2, gas_B2, gas_C2, gas_D2] if Rw is not None else [gas_A2, gas_C2]
    L = [liq_A1(1.0), liq_C1(1.0)]
    G = [g(1.0) for g in gas]
    rows = []
    rows.append([f[0] for f in L] + [-g[0] for g in G])                                       # u_r continuous
    rows.append([f[1] for f in L] + [-g[1] for g in G])                                       # u_z continuous
    rows.append([mu1 * (f[2] + k * f[0]) for f in L] + [-mu2 * (g[2] + k * g[0]) for g in G])  # tangential stress
    rows.append([f[4] + 2 * mu1 * f[3] - (1 - k * k) * f[0] / s for f in L]
                + [-(g[4] + 2 * mu2 * g[3]) for g in G])                                      # normal stress, eta = u_r/s
    if Rw is not None:
        W = [g(float(Rw)) for g in gas]
        rows.append([0.0, 0.0] + [g[0] for g in W])                                           # u_r = 0 at the wall
        rows.append([0.0, 0.0] + [g[1] for g in W])                                           # u_z = 0 at the wall
    return np.array(rows, dtype=float)

def det(s, k, oh, rho2, mu2, Rw=None):
    M = _rows(s, k, oh, rho2, mu2, Rw)
    M = M / np.max(np.abs(M), axis=0)          # column normalisation
    return np.linalg.det(M)

def inviscid_bound(k, rho2, Rw=None):
    """Two-fluid inviscid rate with the wall (the gas potential satisfies u_r = 0 at Rw); an upper bound on the viscous root."""
    k = float(k)
    if Rw is None:
        g = kv(0, k) / kv(1, k)
    else:
        Rk = k * float(Rw)
        g = (kv(0, k) * iv(1, Rk) + iv(0, k) * kv(1, Rk)) / (kv(1, k) * iv(1, Rk) - iv(1, k) * kv(1, Rk))
    return np.sqrt(k * (1 - k * k) / (iv(0, k) / iv(1, k) + float(rho2) * g))

def growth_rate(k, oh, rho2, mu2, Rw=None, n=400):
    """Largest real root in (0, inviscid confined bound]; bracketed by scanning down from the bound."""
    up = inviscid_bound(k, rho2, Rw) * (1 + 1e-9)
    f = lambda s: det(s, k, oh, rho2, mu2, Rw)
    prev_s, prev_v = up, f(up)
    for i in range(1, n + 1):
        si = up * 10.0 ** (-6.0 * i / n)
        vi = f(si)
        if prev_v * vi <= 0:
            return brentq(f, si, prev_s, xtol=1e-15, rtol=1e-15, maxiter=200)
        prev_s, prev_v = si, vi
    return float('nan')


def eigenfunction(r, k, s, oh, rho2, mu2, Rw=None):
    """Radial profiles (u_r, u_z / i) of the eigenmode at the root s, from the null vector of M, using the same scaled basis as _rows (columns A1, C1, A2, [B2], C2, [D2]). Returns (u_r, u_z/i, smallest/largest singular value of the column-normalised M)."""
    M = _rows(s, k, oh, rho2, mu2, Rw)
    scale = np.max(np.abs(M), axis=0)
    _, sv, vt = np.linalg.svd(M / scale)
    c = vt[-1] / scale
    if Rw is None:
        A1, C1, A2, C2 = c; B2 = D2 = 0.0
    else:
        A1, C1, A2, B2, C2, D2 = c
    nu2 = float(mu2) / float(rho2)
    l1 = np.sqrt(k * k + s / oh); l2 = np.sqrt(k * k + s / nu2)
    Rref = 1.0 if Rw is None else float(Rw)
    I = lambda n, a, x, cc: ive(n, a * x) * np.exp(a * x - cc)
    K = lambda n, a, x, cc: kve(n, a * x) * np.exp(cc - a * x)
    r = np.asarray(r, dtype=float); ur = np.empty_like(r); uz = np.empty_like(r)
    i = r <= 1.0; o = ~i; ri, ro = r[i], r[o]
    ur[i] = A1 * k * I(1, k, ri, k) + C1 * I(1, l1, ri, l1)
    uz[i] = A1 * k * I(0, k, ri, k) + C1 * (l1 / k) * I(0, l1, ri, l1)
    ur[o] = -A2 * k * K(1, k, ro, k) + B2 * k * I(1, k, ro, k * Rref) + C2 * K(1, l2, ro, l2) + D2 * I(1, l2, ro, l2 * Rref)
    uz[o] = A2 * k * K(0, k, ro, k) + B2 * k * I(0, k, ro, k * Rref) - C2 * (l2 / k) * K(0, l2, ro, l2) + D2 * (l2 / k) * I(0, l2, ro, l2 * Rref)
    return ur, uz, sv[-1] / sv[0]
