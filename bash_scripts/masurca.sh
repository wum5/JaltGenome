#!/bin/bash

#PBS -N MaSuRCA
#PBS -l nodes=1:ppn=16,walltime=320:00:00,vmem=400gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load gcc/4.9.2
module load perl/5.16.2


cd /N/dc2/projects/jaltomt/GenomeAssembly/Assembly/MaSuRCA_assembly
./assemble.sh
