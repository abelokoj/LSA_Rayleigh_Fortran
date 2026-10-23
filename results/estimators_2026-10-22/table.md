# Every estimator for every case (22 October 2026)

Source: `results/followup_2026-10-22/analysis/summary.csv` and the record reports in `results/followup_2026-10-22/analysis/records/*.txt`. Machine-readable: `all_estimators.csv` (rates and departures). All departures are from the two-fluid viscous rate (Tomotika 1935) for each run's own parameters; the start-time spread is the range of the interface rate when the start of the fit window is moved (upper end fixed), as a percentage of the rate. Velocity DMD is the twin DMD at its selected rank; the rank-scan median and its interquartile range (IQR, as a percentage of the median) are over truncation ranks 2–12. Blank cells: the quantity does not exist for that case (no twin, so no DMD; or a stable case with no growth-rate reference).

## Rule for which number is quoted (Chapter 4 of the dissertation)

1. **One estimator is quoted for every case: the interface-amplitude fit**, with its start-time spread as the uncertainty. It is the only estimator available for every run, it does not need the twin, and its window is chosen by a fixed procedure (opens when the local rate has settled, closes at a/R₀ = 0.1).
2. **Flags, decided by fixed thresholds, not by the result:** a window of less than half an e-fold (the record report already warns at this level) marks the value as poorly constrained; it is reported, but no quantitative claim rests on it alone. A velocity-DMD rank-scan IQR above 1 % marks the velocity DMD as unreliable for that case.
3. **The two DMD estimates are independent checks, always shown, never substituted.** The level-set DMD runs slightly low; Chapter 7 of the dissertation discusses why.
4. **Archived runs** (ε = 0.02: `disp*`, `res*`, `eps*`, `pinch*`, `domain06`) are listed for completeness and the history of the method, not as results.

Consequence for the near-cut-off cases (x = 0.9 and 0.95): the quoted value is the interface fit (+1.4 to +3.8 %), flagged as a short window, with the level-set DMD (−1.4 to −1.9 %) shown beside it as the check; the two bracket the theory.

## Table

| Case | k R₀ | Window (e-folds) | Interface fit | Start-time spread | Level-set DMD | Velocity DMD | Velocity rank-scan median (IQR) | Flags |
|---|---|---|---|---|---|---|---|---|
| amp08_025 | 0.70 | 1.07 | +0.44 % | 0.47 % |  |  |  |  |
| amp08_050 | 0.70 | 0.68 | +0.41 % | 0.52 % |  |  |  |  |
| amp08_100 | 0.70 | 0.42 | +0.00 % | 0.31 % |  |  |  | short window (< 0.5 e-fold) |
| amp32_025 | 0.70 | 1.07 | -0.08 % | 0.42 % |  |  |  |  |
| amp32_050 | 0.70 | 0.64 | -0.10 % | 0.56 % |  |  |  |  |
| amp32_100 | 0.70 | 0.35 | -0.42 % | 0.40 % |  |  |  | short window (< 0.5 e-fold) |
| conf15_16 | 0.70 | 1.31 | -3.03 % | 0.28 % | -3.35 % | -0.20 % | -0.46 % (2.54 %) | velocity DMD: rank-scan IQR > 1 % |
| conf15_32 | 0.70 | 1.48 | -3.11 % | 0.31 % | -3.34 % | -1.31 % | -2.56 % (1.37 %) | velocity DMD: rank-scan IQR > 1 % |
| conf15k050_16 | 0.50 | 1.64 | -3.67 % | 0.03 % | -3.12 % | -4.96 % | -2.36 % (0.04 %) |  |
| conf15k060_16 | 0.60 | 1.53 | -3.27 % | 0.13 % | -3.08 % | -3.11 % | -2.10 % (0.19 %) |  |
| conf15k080_16 | 0.80 | 0.71 | -2.82 % | 0.55 % | -3.66 % | +2.88 % | +8.64 % (15.81 %) | velocity DMD: rank-scan IQR > 1 % |
| conf20_16 | 0.70 | 1.35 | -0.63 % | 0.34 % | -0.91 % | +2.15 % | -0.93 % (3.12 %) | velocity DMD: rank-scan IQR > 1 % |
| conf30_16 | 0.70 | 1.35 | -0.14 % | 0.27 % | -0.45 % | +2.18 % | -0.54 % (2.90 %) | velocity DMD: rank-scan IQR > 1 % |
| disp020 | 0.20 | 1.04 | -6.20 % |  |  |  |  | archived (eps = 0.02; history only) |
| disp030 | 0.30 | 1.04 | -6.18 % |  |  |  |  | archived (eps = 0.02; history only) |
| disp040 | 0.40 | 0.20 | -1.62 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| disp050 | 0.50 | 0.20 | -1.81 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| disp060 | 0.60 | 0.24 | -1.89 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| disp080 | 0.80 | 0.24 | -1.26 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| disp090 | 0.90 | 0.17 | +0.16 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| disp110 | 1.10 | -2.91 |  | 233.38 % |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold); stable: no growth-rate reference |
| dom06_16 | 0.70 | 1.33 | -0.07 % | 0.29 % |  |  |  |  |
| domain06 | 0.70 | 0.24 | -1.77 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| dts125_16 | 0.70 | 1.34 | -0.08 % | 0.29 % |  |  |  |  |
| eps005 | 0.70 | 0.43 | -26.23 % | 11.12 % |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| eps010 | 0.70 | 0.43 | -26.28 % | 11.12 % |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| eps040 | 0.70 | 0.43 | -25.94 % | 11.20 % |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| lin050_16 | 0.70 | 0.72 | -0.14 % | 0.49 % |  |  |  |  |
| lin100_16 | 0.70 | 0.37 | -0.46 % | 0.35 % |  |  |  | short window (< 0.5 e-fold) |
| meas08 | 0.70 | 0.57 | +0.37 % | 0.92 % | -0.80 % | +0.67 % | +0.60 % (0.57 %) |  |
| meas16 | 0.70 | 0.73 | -0.24 % | 0.74 % | -0.66 % | -1.54 % | -1.48 % (9.89 %) | velocity DMD: rank-scan IQR > 1 % |
| meas32 | 0.70 | 0.65 | -0.23 % | 0.82 % | -0.74 % | +2.23 % | +0.84 % (1.27 %) | velocity DMD: rank-scan IQR > 1 % |
| mug10_16 | 0.70 | 1.41 | -0.08 % | 0.26 % |  |  |  |  |
| pinch16 | 0.70 | 0.22 | -2.13 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| pinch32 | 0.70 | 0.22 | -2.00 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| res08 | 0.70 | 0.24 | -1.78 % |  | +0.03 % | -5.65 % | -5.65 % (4.01 %) | archived (eps = 0.02; history only); short window (< 0.5 e-fold); velocity DMD: rank-scan IQR > 1 % |
| res16 | 0.70 | 0.20 | -2.05 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| res32 | 0.70 | 0.20 | -1.95 % |  |  |  |  | archived (eps = 0.02; history only); short window (< 0.5 e-fold) |
| rhog10_16 | 0.70 | 1.33 | -0.11 % | 0.28 % |  |  |  |  |
| stab110_16 | 1.10 | 0.01 |  |  |  |  |  | short window (< 0.5 e-fold); stable: no growth-rate reference |
| sw020_16 | 0.20 | 1.42 | -0.17 % | 0.28 % | -1.31 % | +0.54 % | +0.43 % (0.63 %) |  |
| sw030_16 | 0.30 | 1.58 | -0.22 % | 0.21 % | +0.17 % | +1.28 % | +0.23 % (0.35 %) |  |
| sw040_16 | 0.40 | 1.60 | -0.24 % | 0.20 % | -0.29 % | +1.47 % | +0.29 % (0.26 %) |  |
| sw050_16 | 0.50 | 1.56 | -0.23 % | 0.18 % | -0.29 % | +1.53 % | +0.32 % (0.77 %) |  |
| sw060_16 | 0.60 | 1.48 | -0.17 % | 0.22 % | +0.28 % | +1.75 % | +0.33 % (0.78 %) |  |
| sw070_16 | 0.70 | 1.33 | -0.09 % | 0.29 % | -0.47 % | +3.16 % | +0.35 % (0.74 %) |  |
| sw080_16 | 0.80 | 0.77 | +0.21 % | 0.41 % | -0.70 % | +9.17 % | +10.10 % (2.68 %) | velocity DMD: rank-scan IQR > 1 % |
| sw090_16 | 0.90 | 0.30 | +1.41 % | 0.78 % | -1.39 % | +7.93 % | +8.19 % (0.69 %) | short window (< 0.5 e-fold) |
| sw090e1_16 | 0.90 | 0.24 | +1.54 % | 0.57 % | -1.58 % | +29.93 % | +27.71 % (0.91 %) | short window (< 0.5 e-fold) |
| sw095e1_16 | 0.95 | 0.14 | +3.84 % | 1.39 % | -1.85 % | +16.40 % | +16.40 % (0.00 %) | short window (< 0.5 e-fold) |
| visc010_16 | 0.70 | 1.34 | -0.35 % | 0.40 % | -0.73 % | -0.17 % | -0.09 % (0.13 %) |  |
| visc030_16 | 0.70 | 1.33 | -0.88 % | 0.37 % | -1.24 % | +0.23 % | -0.39 % (0.61 %) |  |
| visc100_16 | 0.70 | 1.50 | -2.03 % | 0.28 % | -2.30 % | -1.67 % | -2.30 % (0.30 %) |  |
| visc100r6_16 | 0.70 | 1.48 | -0.33 % | 0.24 % | -0.61 % | +0.69 % | -0.15 % (0.27 %) |  |
| visc100r8_16 | 0.70 | 1.51 | -0.11 % | 0.23 % | -0.37 % | +0.23 % | -0.05 % (0.18 %) |  |
