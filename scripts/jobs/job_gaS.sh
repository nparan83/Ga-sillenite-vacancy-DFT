#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_S40
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
cd /arf/scratch/nparan/SilenitB/02_O40/gaS_rel || exit 1
mpirun pw.x -npool 2 -in run.in > run.out || { grep -q '^!' run.out || mpirun pw.x -npool 1 -in run.in > run.out; }
