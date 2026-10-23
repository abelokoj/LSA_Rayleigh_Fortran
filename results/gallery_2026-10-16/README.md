# Figure gallery of the October 2026 campaign (2026-10-16)

This is the October counterpart of the September `all_figures.pdf` and GIFs (8 cells per radius, ε = 0.02, archived with the September analysis package in the solver fork `abelokoj/amrex_implicit_interfaces`, commit `dddb827c`, `NS_fluids_lib/LSA_Rayleigh_Fortran/figures`). It was drawn by `plotting/plot_campaign_gallery.py` from the snapshot now in `results/campaign_2026-10-09`, the two-fluid records of `2026-10-16_two_fluid`, and the velocity export of meas16. No plotfile was needed. Every static figure is in PNG, PDF and SVG; `all_figures.pdf` has one per page.

| September figure | October figure here | Source data |
|---|---|---|
| `growth_kr0_0p7` | `gallery_growth_sw070_16`, with the final fit window, Chandrasekhar and two-fluid viscous slopes | `records/sw070_16.dat`, record CSV |
| `local_sigma_kr0_0p7` | `gallery_local_rate_sw070_16` | same |
| `interface_kr0_0p7` | `gallery_interface_sw070_16` | `bundles/sw070_16/profiles.dat` |
| `fig3_proportionality` | `gallery_linearity` (ε = 2.5e-3, 5e-3, 1e-2 at 16 cells) | records sw070_16, lin050_16, lin100_16 |
| `fig4_dispersion`, `fig5_dispersion_measured` | `campaign_dispersion`; new `gallery_sweep_growth` | `analysis/summary.csv`, sweep records |
| `field_y_velocity_kr0_0p7` | `gallery_field_meas16`: the perturbation (twin difference) rather than the full field, each panel scaled by its maximum | velocity export of meas16 |
| `pinchoff_sequence`, `_extrema`, `_extrapolation` | same names, from pinch16 (16 cells) | `bundles/pinch16` |
| `ohnesorge_contour`, `_surface`, `_peak_drift` | same names (theory only) | none |
| `schematic_twin_run` | same | none |
| `anim_interface.gif` | `anim_interface_sw070_16.gif` (linear range) and `anim_interface_pinch16.gif` (through break-up), each with an `.html` slider | profile bundles |
| `anim_field_y_velocity.gif` | `anim_field_meas16.gif` (perturbation axial velocity with the interface) | velocity export of meas16 |
| `fig1_validation`, `fig2_replication` | not redrawn: they compared solver builds of September, which the campaign does not repeat | — |
| (none) | `campaign_stable_band` | stab110_16 |

Not yet possible from workstation data:

- **Full field of the perturbed run at 16 or 32 cells.** The gallery shows the twin difference instead. The full field needs the plotfiles on the cluster, or `dump_profiles` field slices written into the bundles.
- **Break-up at 32 cells.** Run `plot_pinchoff.py --bundle .../bundles/pinch32 --cells-per-r0 32` for the mesh comparison of the break-up time.

To redraw (from the top of the repository, about 6 minutes, most of it in the GIFs; add `--no-gif` to skip them):

```
python3 plotting/plot_campaign_gallery.py results/campaign_2026-10-09 --two-fluid results/two_fluid_2026-10-16/records --export results/exports/meas16 --export-case meas16 --outdir figures/gallery
```

`results/exports/meas16/` is the lossless velocity export of the twin difference of `meas16` (`perturbation.bin`, 24 MB), included so that `gallery_field_meas16`, the perturbation-field figure of Chapter 5, can be redrawn.
