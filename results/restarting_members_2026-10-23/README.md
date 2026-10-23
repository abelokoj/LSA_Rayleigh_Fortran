# Restarting-driver members, per member (2026-10-23)

These files give the per-member amplitude and growth of the solver's restarting Krylov driver at k r0 = 0.7 and 16 cells, for a figure of the restarting variant. They were computed from the lossless exports of the members (`krylov_test/analysis/<label>/export_krylov/`, written by `krylov_test/krylov_run.sh`; 9–13 MB each, not included), which hold the level set and both velocity components of every member. No new cluster run was needed.

The two runs are:

| File | Run | Members | τ |
|---|---|---|---|
| `kry070_16b_members.csv` | second trial, the one used in the thesis cross-check | 40 | 0.2500 |
| `kry070_16_members.csv` | first trial; its first member has the solver's hard-coded 2-cell level-set amplitude, so its early members are nonlinear | 60 | 0.0808 |

The columns are:

- `member`, `t = (member − 1) τ`;
- `ls_amplitude`: the axial Fourier amplitude, at the imposed wavenumber, of the level-set perturbation averaged over the band |r − 1| < 2 Δr. In the band φ' ≈ −η', so this is the interface displacement.
- `velocity_norm`: the volume-weighted L2 norm of (u_r, u_z);
- `rate_ls`, `rate_velocity`: ln(value_k / value_(k−1)) / τ between consecutive members;
- `two_fluid_viscous`: the unconfined two-fluid viscous rate, 0.3287, for reference.

The last-member rates are, for kry070_16b, 0.3248 (level set) and 0.3281 (velocity); for kry070_16, 0.2971 and 0.3270. These member exports come from the older driver, which did not rescale its members, so the ratios between members are growth factors.
