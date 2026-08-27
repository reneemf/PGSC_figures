# Fig2 — PGS simulations

Self-contained simulations comparing **PGS**, **ampPGS**, and **PGSC** across genetic-architecture parameters. Each sweep varies one parameter and records R².

## Requirements

R 4.3.1 with `dplyr`

## Run

```bash
export PGSC_HOME="/path/to/Figures_repo/"

sbatch run_sims.sh          # 1. SLURM array (--array=1-1000); writes Rdata/<sweep>/*_iter_*.Rdata
Rscript collapse_results.R  # 2. collapse per-iteration files -> Rdata/<sweep>/*_collapsed.Rdata
Rscript panel_summary_plots.R  # 3. render panels -> figs/*.png
```

Edit the `#SBATCH -o/-e` log paths in `run_sims.sh` for your environment.

## Simulation parameters

| Script | Parameter varied |
|---|---|
| `h2.R` | Additive heritability (h²) |
| `h2_gxe.R` | Uncoordinated GxC heritability |
| `h2_coord.R` | Coordinated GxC heritability |
| `eta.R` | Heteroskedasticity (η) |
| `poly.R` | Polygenicity (S_caus) |
| `pop.R` | Sample size (N) |
| `lam_box.R` | Box-Cox scaling (λ) |
| `alpha.R` | Context main effect (α) |
| `c_prob.R` | Context distribution (P(C=1)) |

Each script sources `defaults.R` (parameters), `functions.R` (simulation + PGS builders), and `build_sim.R`.
Most scripts loop over two GxC backgrounds, `gxe0` (none) and `gxenz` (0.05). `panel_summary_plots.R` renders panels A–D to `figs/`:

| Figure(s) | Contents |
|---|---|
| `panel_A.png`, `panel_A_logeta.png` | R² %Change vs GxC heritability, proportion coordinated, and heteroskedasticity (zero GxC background; `_logeta` uses a log η x-axis) |
| `panel_B.png` | Same three sweeps as absolute R² |
| `panel_C_diff.png`, `panel_C_R2.png` | 2×4 grid — additive h², polygenicity, sample size, context imbalance — with the two GxC backgrounds as rows (R² %Change and R²) |
| `panel_D_diff.png`, `panel_D_R2.png` | 2×4 grid — Box-Cox scaling and main effect of C (log / true / squared scales) — with the two GxC backgrounds as rows (R² %Change and R²) |

