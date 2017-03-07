#!/bin/bash

#PBS  -o BATCH_OUTPUT 
#PBS  -e BATCH_ERRORS 

#PBS -lnodes=1:ppn=16:gpu,walltime=4:00:00

# Change to directory from which qsub was executed 
cd $PBS_O_WORKDIR


#Note: Reserve full nodes when using the gpu(s)
module load R
R CMD BATCH --vanilla paschold.R 