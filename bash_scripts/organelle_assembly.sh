#!/bin/bash

#PBS -N Organelle_MT
#PBS -l nodes=1:ppn=8,walltime=120:00:00,vmem=320gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load perl/5.16.2
module load samtools
module load gcc/4.9.2
module load bedtools/2.26.0 
module load celera/8.3rc2
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
module load bioperl/1.6.1


PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/Organelle_PBA/seqtk
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/Organelle_PBA
export BLASR_PATH=/N/dc2/projects/jaltomt/Softwares/Organelle_PBA/SSPACE-LongRead_v1-1
export SPRAI_PATH=/N/dc2/projects/jaltomt/Softwares/Organelle_PBA/sprai-0.9.9.23/bin
export SSPACELONG_PATH=/N/dc2/projects/jaltomt/Softwares/Organelle_PBA/SSPACE-LongRead_v1-1
cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
cpanm Perl4::CoreLibs
cpanm IPC/Run.pm


cd /N/dc2/projects/jaltomt/GenomeAssembly/Assembly/Organelle_Assembly

# Chloroplast
OrganelleRef_PBA -i DHS1_filtered_subreads.fastq -r SLcp.fa -o chloroplast -b '-nproc=8' -s 'num_threads=8'

# Mitochondria
OrganelleRef_PBA -i DHS1_filtered_subreads.fastq -r NTmt.fa -o mitochondira -b '-nproc=8' -s 'num_threads=8'

