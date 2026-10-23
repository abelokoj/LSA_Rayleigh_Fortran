"""
check_inviscid_limit.py
=======================
Independent check of the wall rows of confined_dispersion.py (Appendix D of the dissertation).

As both viscosities tend to zero (fixed viscosity ratio), the viscous confined growth rate must tend to the closed-form inviscid confined rate

    s^2 = k (1 - k^2) / ( I0(k)/I1(k) + rho2 * G ),   G = [K0(k) I1(k Rw) + I0(k) K1(k Rw)] / [K1(k) I1(k Rw) - I1(k) K1(k Rw)],

derived separately from the potential-flow problem with u_r = 0 at r = Rw (capillary units). The test uses an ambient as dense as the liquid (rho2 = 1), so that the wall changes the inviscid rate substantially and the limit discriminates between the confined and the unconfined formula; the algebra being tested does not depend on the parameter values. The departure must shrink steadily (boundary layers of thickness ~ sqrt(nu)) towards the confined value, not the unconfined one.
"""
import numpy as np
from confined_dispersion import growth_rate, inviscid_bound

rho2, ratio = 1.0, 0.5
print(f"ambient density {rho2}, viscosity ratio mu2/mu1 = {ratio}")
for k in (0.5, 0.7, 0.9):
    for Rw in (1.5, 2.0, 4.0):
        sc, su = inviscid_bound(k, rho2, Rw), inviscid_bound(k, rho2, None)
        line = f"k={k} Rw={Rw}: inviscid confined {sc:.6f} (unconfined {su:.6f}, wall effect {100*(sc/su-1):+.2f} %); viscous - confined inviscid:"
        prev = None
        for oh in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10):
            s = growth_rate(k, oh, rho2, ratio * oh, Rw)
            d = s / sc - 1
            line += f"  oh={oh:.0e}: {d:+.2e}"
        print(line)
