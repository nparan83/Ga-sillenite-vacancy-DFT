#!/bin/bash
#SBATCH -A nparan
#SBATCH -p hamsi
#SBATCH -J gaB_ref
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=04:00:00
#SBATCH -o slurm-%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
cd /arf/scratch/nparan/SilenitB || exit 1
for f in aGa bGa2O3; do
  np=$([ $f = aGa ] && echo 8 || echo 4)
  ( cd 06_ref/${f}_vc && [ ! -s run.out ] && mpirun pw.x -npool $np -in run.in > run.out )
  grep -q 'JOB DONE' 06_ref/${f}_vc/run.out && grep -q 'bfgs converged' 06_ref/${f}_vc/run.out || { echo "### $f vc-relax yakinsamadi -> scf yok"; continue; }
  python3 99_analiz/ga_ref.py scf $f || continue
  ( cd 06_ref/${f}_scf && mpirun pw.x -npool $np -in run.in > run.out )
  if grep -q 'JOB DONE' 06_ref/${f}_scf/run.out; then rm -f 06_ref/${f}_vc/out/*.wfc[0-9]* 06_ref/${f}_vc/out/*.save/wfc*.dat 06_ref/${f}_scf/out/*.wfc[0-9]* 06_ref/${f}_scf/out/*.save/wfc*.dat; echo "temizlik $f"; fi
done
