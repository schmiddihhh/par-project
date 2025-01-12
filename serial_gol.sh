#!/bin/bash
#SBATCH --job-name=job_task_1
#SBATCH --output=task_1.out
#SBATCH --error=task_1.err
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=1
#SBATCH --time=10:00

# Load any necessary modules (if needed)
# module load module_name

# Enter your executable commands here
# Execute the compiled program
module load compiler/gcc/
module load mpi/openmpi/
mpicxx task_1.cpp io.cpp
mpirun -n 1 a.out
