#!/bin/bash
#SBATCH --job-name=SPECTRAL_trying_k
#SBATCH --nodes=1
#SBATCH --ntasks=128   #128   # EITHER 128, OR 1 OR 1, OR 81 DEPENDING ON FILE
#SBATCH --cpus-per-task=1           # 1 CPU per process
#SBATCH --time=2-00:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ykv210@nyu.edu
#SBATCH --mem=400GB


module load matlab

matlab -nodisplay -nosplash -r "run('run_casa.m'); exit"

