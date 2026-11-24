#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_eos
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=12:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
cd /arf/scratch/nparan/SilenitB/04_eos/$1 || exit 1
[ -s run.out ] && { echo "### run.out zaten var -> DUR"; exit 1; }
mpirun pw.x -npool 2 -in run.in > run.out
if grep -q 'JOB DONE' run.out && { ! grep -qi "calculation *= *'relax'" run.in || grep -q 'bfgs converged' run.out; }; then
  rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik: gecici dalga fonksiyonlari silindi"
fi
