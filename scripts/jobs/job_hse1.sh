#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_hse1
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=16:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
cd /arf/scratch/nparan/SilenitB/07_hse/ga${1}_hse || exit 1
[ -s run.out ] && { echo "### run.out zaten var -> DUR"; exit 1; }
mpirun pw.x -npool 1 -in run.in > run.out
grep -q 'JOB DONE' run.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik hse $1"
