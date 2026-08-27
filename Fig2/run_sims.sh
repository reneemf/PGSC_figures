#!/usr/bin/env bash

#SBATCH -J run_sims
#SBATCH --mem=80GB
#SBATCH --time=10:00:00
#SBATCH --partition=tier1q
#SBATCH -o /home/<username>/slurm_outputs/run_sims/run_sims_%a_%A.out
#SBATCH -e /home/<username>/slurm_outputs/run_sims/run_sims_%a_%A.err
#SBATCH --array=1-1000

start_time=$(date +%s)  
echo "Date run: $(date)"

# Load required modules
module load gcc/12.1.0
module load R/4.3.1

# Map array task ID to iteration index
iter_idx="${SLURM_ARRAY_TASK_ID}"

# h2_coord and h2_gxe: no GxC-background split
echo "Running: Rscript h2_coord.R ${iter_idx}"
Rscript h2_coord.R "${iter_idx}"
echo "Running: Rscript h2_gxe.R ${iter_idx}"
Rscript h2_gxe.R "${iter_idx}"

# h2, pop, poly, lam_box, eta: w GxC-background split
echo "Running: Rscript h2.R ${iter_idx}"
Rscript h2.R "${iter_idx}"
echo "Running: Rscript pop.R ${iter_idx}"
Rscript pop.R "${iter_idx}"
echo "Running: Rscript poly.R ${iter_idx}"
Rscript poly.R "${iter_idx}"
echo "Running: Rscript lam_box.R ${iter_idx}"
Rscript lam_box.R "${iter_idx}"
echo "Running: Rscript alpha.R ${iter_idx}"
Rscript alpha.R "${iter_idx}"
echo "Running: Rscript eta.R ${iter_idx}"
Rscript eta.R "${iter_idx}"
echo "Running: Rscript c_prob.R ${iter_idx}"
Rscript c_prob.R "${iter_idx}"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "Elapsed time: ${elapsed} seconds"


