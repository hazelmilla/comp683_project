#! /bin/bash -login
#SBATCH -J impute
#SBATCH -t 2-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -p general
#SBATCH --mem=128gb
#SBATCH --mail-user=dpguilba@email.unc.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -o impute_missing_noval.out

set -e

# Activate conda
source /nas/longleaf/home/dpguilba/miniforge3/etc/profile.d/conda.sh

conda create -n impute python=3.11 -y

conda activate impute

mamba install -c conda-forge pandas numpy scipy statsmodels scikit-learn matplotlib pyarrow -y

# Optional: ensure libraries use all cores
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8
export NUMEXPR_NUM_THREADS=8

python -u impute_missing_noval.py
