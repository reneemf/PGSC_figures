# Fig2 — PGS simulations

Self-contained simulations comparing **PGS**, **ampPGS**, and **PGSC** across genetic-architecture parameters. Each sweep varies one parameter and records R²; no external data needed.

## Requirements

R 4.3.1 with `dplyr`, GCC 12.1.0, and SLURM (for `run_sims.sh`).

## Run

```bash
export PGSC_HOME="/path/to/Figures_repo/"

sbatch run_sims.sh          # 1. SLURM array (--array=1-1000); writes Rdata/<sweep>/*_iter_*.Rdata
Rscript collapse_results.R  # 2. collapse per-iteration files -> Rdata/<sweep>/*_collapsed.Rdata
Rscript panel_summary_plots.R  # 3. render panels -> figs/*.png
```

Edit the `#SBATCH -o/-e` log paths in `run_sims.sh` for your environment. As shipped, `run_sims.sh` runs `alpha.R` and `eta.R`; uncomment the other sweep blocks to enable them.

## Sweeps

Each sweep sources `defaults.R` (parameters), `functions.R` (simulation + PGS builders), and `build_sim.R` (`run_sim()`).

| Script | Parameter varied |
|---|---|
| `h2.R` | Additive heritability (h²) |
| `h2_gxe.R` | Uncoordinated GxC heritability |
| `h2_coord.R` | Coordinated GxC proportion |
| `eta.R` | Heteroskedasticity (η) |
| `poly.R` | Polygenicity (S_caus) |
| `pop.R` | Sample size (N) |
| `lam_box.R` | Box-Cox scaling (λ) |
| `alpha.R` | Context main effect (α) |
| `c_prob.R` | Context prevalence (P(C=1)) |

Most sweeps loop over two GxC backgrounds, `gxe0` (none) and `gxenz` (nonzero); `panel_summary_plots.R` writes `panel_{A,B,C,D}_{gxe0,gxenz}.png`. Key defaults (`defaults.R`): `N=1e4`, `S=1e4`, `S_caus=1e3`, `h2_add=0.3`, `h2_gxe=0.2`, 1000 iterations.
