#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_yak
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=06:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
R=/arf/scratch/nparan/SilenitB/01_yakinsama
for n in S_e080k2 S_e070k2 S_e090k2 S_e100k2 S_e080k1 S_e080k3 B_e080k2 A_e080k2 B_e100k2 A_e100k2; do
  [ -s $R/$n/run.in ] || continue
  np=2; case $n in *k1) np=1;; esac
  ( cd $R/$n && { mpirun pw.x -npool $np -in run.in > run.out || mpirun pw.x -npool 1 -in run.in > run.out; } \
    && grep -q 'JOB DONE' run.out && rm -f out/*.save/wfc* )
done
