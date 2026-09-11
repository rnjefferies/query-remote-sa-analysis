# Query-Remote-SA Analysis (Study 1)

Reproducible analysis package for **Chapters 3 and 4** of the thesis *Towards Passive
Measurement of Situation Awareness in the Remote Operation of Live Vehicles*. Study 1 had
24 operators remotely drive a live electric vehicle in a within-subjects design across three
situation-awareness query methods: **No Query (A)**, **Realtime (B, answering while driving)**,
and **Offline (C, stopping to answer)**.

- **Chapter 3 (offline versus real-time SA):** subjective workload (NASA-TLX), objective and
  subjective situation awareness (SART), and remote-driving performance (steering, throttle,
  braking, collisions, speed).
- **Chapter 4 (confidence):** confidence calibration, overconfident errors, and the type-2
  (metacognitive) discrimination of correct from incorrect answers.

This is the **analysis** half of Study 1. The data-collection and processing tools are in the
companion repository
**[eyetrack-remote-sa-tools](https://github.com/rnjefferies/eyetrack-remote-sa-tools)**; the
Study 2 analyses are in **eyetrack-remote-sa-analysis**.

> **Data included.** Only de-identified, analysis-ready CSVs are distributed (pseudonymous
> operator codes such as `01_RT`; no raw survey exports, consent, signatures, names, emails,
> IP addresses, or GPS). The raw Qualtrics export and consent records are **not** distributed.
> Per-participant demographics are also withheld (the analyses do not use them); the aggregate
> demographic counts are provided in `data/demographics_summary.csv`.

## Repository layout

```
ch3_offline_vs_realtime/    workload, SART, objective SA, and driving-performance analyses
ch4_confidence/             confidence calibration, overconfidence, type-2 AUROC
data/                       de-identified input CSVs
outputs/                    (ships empty) result tables written when you run the scripts
renv.lock                   pinned R package versions (R 4.2.1)
```

## Quick start

Run from the **repository root** in R (developed on R 4.2.1):

```r
install.packages("renv")     # if not already installed
renv::restore()              # installs the exact pinned package versions from renv.lock
```

```bash
Rscript ch3_offline_vs_realtime/ch3_offline_vs_realtime.R
Rscript ch4_confidence/confidence_pipeline.R      # writes result tables to outputs/
```

- Run every command from the repository root so the `data/` and `outputs/` paths resolve.
- The Ch3 script prints each computed statistic to the console as it runs.
- All analyses are deterministic.

## Chapter 3: offline versus real-time SA

`ch3_offline_vs_realtime.R` reproduces, from the de-identified CSVs:

| Analysis | Reported result |
|---|---|
| Overall workload (NASA-TLX), RM-ANOVA | F = 15.29, partial eta-squared = .17 |
| Subscale Friedmans (Mental / Frustration / Effort) | chi-square = 24.30 / 15.50 / 11.60 |
| SART Demand / Total, RM-ANOVA | F = 14.86 / 3.82 |
| Objective SA, paired t-test (Realtime vs Offline) | t(23) = 5.76, mean difference 12.81 |
| Throttle-reversal rate, RM-ANOVA | F = 34.37 |
| Fine / Major steering-reversal rate, Friedman | chi-square = 18.6 / 32.2 |
| Braking rate (per minute), Friedman | chi-square = 22.9 |
| Collision rate (per minute), Friedman | chi-square = 0.64 (ns) |
| Normalised speed mean / SD, RM-ANOVA | F = 1.16 (ns) / 8.34 |

## Chapter 4: confidence

`confidence_pipeline.R` is a CSV-in, CSV-out pipeline that rebuilds the query-level and
operator-level confidence measures, calibration and resolution indices (with cluster bootstrap
CIs), overconfident-error counts, and the type-2 AUROC comparison (real-time vs offline). All
result tables are written to `outputs/`.

## Citation

If you use this code or data, please cite the thesis (R. Jefferies, *Towards Passive Measurement
of Situation Awareness in the Remote Operation of Live Vehicles*). Full citation and DOI to be
added on deposit.

## Licence

Code released under the MIT Licence (see [`LICENSE`](LICENSE)). The de-identified data are shared
for reproduction of the reported analyses.
