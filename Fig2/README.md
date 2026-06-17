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
Most scripts loop over two GxC backgrounds, `gxe0` (none) and `gxenz` (nonzero); `panel_summary_plots.R` writes `panel_{A,B,C,D}_{gxe0,gxenz}.png`. 
