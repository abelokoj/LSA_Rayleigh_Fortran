# Archived amplitude records

These are the interface-amplitude records of the September 2026 study, copied on 2026-10-07 from the September analysis package (`NS_fluids_lib/LSA_Rayleigh_Fortran/results/records/` of the solver fork `abelokoj/amrex_implicit_interfaces`, commit `dddb827c`). Each has three columns: time, mode amplitude and mean radius. They predate the six-column format written by `bin/dump_amplitude` since 2026-10-07.

The records are kept here so that every interface growth rate of the archived study can be recomputed without the cluster. `window_analysis.txt` is the output of `bin/run_lsa record` on each of them, with the defaults (a/r0 ≤ 0.1, sliding width 0.8, settling tolerance 0.5 %). Its columns are label, k r0, window start and end, samples, growth rate, e-folds spanned, settled (T/F), bound reached (T/F), exact viscous rate and percentage departure from it.

What it shows:

- **The small-amplitude runs:** only `meas08`, `meas16` and `meas32` (ε = 5 × 10⁻³) have a settled window of more than half an e-fold before a/r0 = 0.1. Their rates are 0.3299, 0.3279 and 0.3279, within 0.5 % of the exact viscous value 0.32937.
- **The ε = 2 × 10⁻² runs:** every one (`res*`, `disp040`–`disp090`, `domain06`, `pinch*`) reaches the amplitude bound after at most 0.24 e-folds of settled growth. Their rates, about 2 % low, are finite-window values. `disp020` and `disp030` never settle before the bound.
- **The amplitude sweep:** `eps005`, `eps010` and `eps040` stop at t = 3.6, inside the startup transient. They support the proportionality test only at early times, and they give no growth rate.
- **The stable control:** `disp110` cannot be fitted with an exponential, because its magnitude folds a damped oscillation. The new campaign's `stab110_16` writes the signed coefficients needed for `plotting/fit_stable_band.py`.

These results decide what the October 2026 campaign (`campaign_*.txt`) runs and what it does not repeat.
