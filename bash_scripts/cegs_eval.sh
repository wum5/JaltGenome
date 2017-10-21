#!/bin/bash

#PBS -N CEGs_EVAL
#PBS -l nodes=1:ppn=2,walltime=96:00:00,vmem=32gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


## for BUSCO
module load python
module load gcc/4.9.2
module load hmmer/3.1b2 
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/augustus-3.2.3/scripts/
export AUGUSTUS_CONFIG_PATH=/N/dc2/projects/jaltomt/Softwares/augustus-3.2.3/config
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/augustus-3.2.3/bin
OD=/N/dc2/projects/jaltomt/Softwares/busco
cd /N/dc2/projects/jaltomt/GenomeAssembly/Assembly/MaSuRCA_assembly

# run against eukaryote CEGs                                               
python $OD/BUSCO.py -i final_assembly.fasta -o eukaryo_CEGs -l $OD/eukaryota_odb9/ -m genome -c 4 -f
 
# run against plant-specific CEGs 
python $OD/BUSCO.py -i final_assembly.fasta -o plant_CEGs -l $OD/embryophyta_odb9/ -m genome -c 8 -f


## for CEGMA 
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/CEGMA_v2/geneid/bin
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/CEGMA_v2/wise2.4.1/src/bin
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/CEGMA_v2/bin
export WISECONFIGDIR=/N/dc2/projects/jaltomt/Softwares/CEGMA_v2/wise2.4.1/wisecfg
export PERL5LIB=/N/dc2/projects/jaltomt/Softwares/CEGMA_v2/lib
export CEGMA=/N/dc2/projects/jaltomt/Softwares/CEGMA_v2
cd /N/dc2/projects/jaltomt/GenomeAssembly/Assembly/MaSuRCA_assembly/CEGMA_Eukaryo
  
cegma -g ../final_assembly.fasta -o jaltomata -T 2

