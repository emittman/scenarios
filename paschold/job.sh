#!/bin/bash

salloc -N 1 -n 16 -p gpu -t 4:00:00

# Change to directory from which qsub was executed 
cd $PBS_O_WORKDIR
module load R
srun R --vanilla CMD BATCH paschold.R
