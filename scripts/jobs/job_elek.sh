#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_elek
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=06:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
y=$1; E=/arf/scratch/nparan/SilenitB/05_elek
cd $E/ga${y}_nscf || exit 1
[ -s run.out ] || mpirun pw.x -npool 4 -in run.in > run.out
grep -q 'JOB DONE' run.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik nscf $y"
cd $E/ga${y}_paw || exit 1
[ -s run.out ] || mpirun pw.x -npool 2 -in run.in > run.out
grep -q 'JOB DONE' run.out || { echo "### PAW scf $y bitmedi -> projwfc yok"; exit 1; }
mkdir -p pdos; mpirun projwfc.x -npool 2 -in proj.in > proj.out
grep -q 'JOB DONE' proj.out && rm -f out/*.wfc[0-9]* out/*.save/wfc* && echo "temizlik paw $y"
