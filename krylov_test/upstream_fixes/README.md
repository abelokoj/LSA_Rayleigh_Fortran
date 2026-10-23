# Check of the proposed fixes to the solver's `fort_process_LSA`

`test_fort_process_LSA_fixes.f90` repeats, on one Gram matrix, the steps of `fort_process_LSA` (msussman42/amrex_implicit_interfaces, `PROB_CPP_PARMS.F90`, commit `3752065`):

- the SVD solve of eq. 67;
- the companion matrix;
- the eigenvalues by `ZGGEV`.

It then selects the mode two ways:

- **original:** every singular value above zero is kept, and the eigenvalue with the largest real part is chosen;
- **corrected:** singular values below `rcond` × the largest are dropped, and the mode chosen is the one of largest Ritz amplitude among those with ln|λ| > δ and |arg λ| ≤ 0.1 ln|λ|.

The Gram matrices are:

- `gram_synthetic.txt`: a synthetic chain of 20 members with one growing mode, s = 0.3287;
- `g_band_*.txt`: members of the restarting-driver trial `kry070_16b` (k r0 = 0.7, τ = 0.25, level set and velocity). They are weighted as upstream does, with 1 within 6 Δx of the interface and 1e-4 outside. The suffix is the first member used: 20 means members 21–40, 30 means 31–40, 10 means 11–30.
- `g_rweight_20.txt`: members 21–40 of `kry070_16b`, all fields, radius weight.

Build and run (system LAPACK):

```
gfortran -O2 test_fort_process_LSA_fixes.f90 -o t -llapack -lblas
./t g_band_20.txt 20 1e-6 1e-2      # file, number of members, rcond, delta
```

Results, as s = ln|λ|/τ (exact viscous 0.3294, two-fluid viscous 0.3287):

| Gram | original | corrected, rcond 1e-10, δ = 0 | corrected, rcond 1e-6, δ = 1e-2 |
|---|---|---|---|
| synthetic | 0.3287 | 0.3287 | 0.3287 |
| band, members 21–40 | 0.7221 (spurious) | 0.3204 | 0.3280 |
| band, members 31–40 | 0.3417 | 0.0023 (neutral mode) | 0.3285 |
| band, members 11–30 | 0.3251 | 0.3251 | 0.3256 |
| radius weight, members 21–40 | 0.6283 (spurious) | 0.3111 | 0.3275 |

The cutoff and the amplitude rule are both needed. The amplitude rule alone, with δ = 0, can pick a near-neutral mode (|λ| ≈ 1.0006) whose amplitude is large but which does not grow. The neutral floor δ removes it.

These chains come from the older driver, which did not rescale its members. They therefore test the solve and the selection only, not the rescaling.
