#!/bin/bash


module load R
srun R --vanilla CMD BATCH paschold.R
