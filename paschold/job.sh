#!/bin/bash
salloc -N 1 -n 16 -p gpu -t 4:00:00
module load R
ls
R --vanilla CMD BATCH paschold.R
