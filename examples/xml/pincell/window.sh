#!/bin/bash

#SBATCH -N 1
#SBATCH -n 8            # 1 cores
#SBATCH -J godiva      # sensible name for the job
#SBATCH -p sched_mit_nse # NSE partition
#SBATCH -o godiva_1group_%A.out
#SBATCH --mem 100000

module use /home/software/nse/modulefiles/ &&  module load blas/OpenBlas/gnu-6.2.0/1 &&  module load openmc/gnu-6.2.0/0.9.0 &&  module load cmake/gnu-6.2.0/3.8.1
module load python/3.6.3

~/code/gptmc/openmc/build/bin/openmc
