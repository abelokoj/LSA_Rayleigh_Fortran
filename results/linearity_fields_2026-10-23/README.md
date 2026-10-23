# Field-direction linearity test (2026-10-23)

At k R₀ = 0.7, the twin-difference fields of three amplitudes, ε = 2.5e-3, 5e-3 and 1e-2, were compared after scaling each by its ε. The comparison was made at matched times, up to the time at which the largest amplitude reaches a/R₀ = 0.1. Three measures were used, each against the smallest amplitude:

- the relative difference of the scaled fields;
- the angle between them;
- the ratio of their norms.

The cases were run on three meshes, each with its own twin: 8 cells (`amp08_*`), 16 cells (`sw070_16`, `lin050_16`, `lin100_16`) and 32 cells (`amp32_*`). The 16-cell runs `lin050_16` and `lin100_16` were re-run on 22–23 October to regenerate their plotfiles. Their amplitude records are identical to the October originals.

Produced on the cluster by `scripts/linearity_fields.sh`, which calls `plotting/linearity_fields.py`. The band-2 results use `--band 2`.

Largest values over the window, given as ε = 5e-3 / ε = 1e-2:

| Mesh | Observable | Relative difference | Angle | Norm ratio |
|---|---|---|---|---|
| 8 | velocity | 3.96 / 8.52 % | 2.23 / 4.84° | 0.88 / 1.80 % |
| 8 | level set, band 2Δx | 16.04 / 17.01 % | 9.11 / 9.65° | 1.28 / 1.44 % |
| 8 | level set, whole domain | 51.6 / 75.3 % | 22.3 / 32.1° | 49.2 / 73.3 % |
| 16 | velocity | 4.13 / 8.60 % | 2.33 / 4.90° | 0.85 / 1.67 % |
| 16 | level set, band 2Δx | 5.73 / 6.10 % | 3.28 / 3.49° | 0.16 / 0.19 % |
| 16 | level set, whole domain | 58.3 / 77.5 % | 32.2 / 46.1° | 48.7 / 73.1 % |
| 32 | velocity | 3.58 / 7.22 % | 2.03 / 4.13° | 0.61 / 1.04 % |
| 32 | level set, band 2Δx | 2.03 / 2.61 % | 1.16 / 1.50° | 0.15 / 0.20 % |
| 32 | level set, whole domain | 57.4 / 74.6 % | 34.8 / 46.2° | 48.6 / 72.7 % |

Reading:

- **Velocity and narrow-band level set: linear in direction as well as in size.** The level-set residual decreases with the mesh (17 → 6 → 2.6 %) and hardly depends on ε, so it is discretization error of the interface reconstruction. The velocity residual doubles with ε and hardly depends on the mesh, which is the quadratic remainder near the end of the window.
- **Whole-domain level set: not proportional to ε, already at t = 0.** The solver stores the level set only within about 0.3 R₀ of the interface and holds it at a fixed value beyond. The cells at the edge of that region jump between the stored and the fixed value by up to 0.14, almost independently of ε. They dominate the whole-domain norm, which gives the norm mismatch 1 − ε_ref/ε = 50 % and 75 %. This is why the level-set observable is restricted to a narrow band.

Files: `lin<mesh>_export.*` (velocity), `lin<mesh>_export_ls.*` (level set, whole domain) and `lin<mesh>_export_ls_band2.*` (level set, band). The `.txt` files hold the per-time tables and the `.csv` files the same tables. The two `.tgz` files are the archives as fetched from the cluster.
