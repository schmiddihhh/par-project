#!/bin/bash
#SBATCH --job-name=job_task_2
#SBATCH --output=task_2.out
#SBATCH --error=task_2.err
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00

# Load any necessary modules (if needed)
# module load module_name

# Enter your executable commands here
# Execute the compiled program
module load compiler/gcc/
module load mpi/openmpi/
mpicxx task_2.cpp io.cpp

mpirun -n 32 a.out --weak_scaling_study 1 100 >> weak_scaling_study.txt
