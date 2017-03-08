#!/bin/bash


module load R
srun -N 1 -n 16 -p gpu -t 4:00:00 R --vanilla CMD BATCH paschold.R
