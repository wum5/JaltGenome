#!/bin/bash

#PBS -N CAFE
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load python
module load biopython

PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/CAFE/release
cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs

## filter input data for CAFE
python ../Scripts/remove_redundancy.py -i output/pairs/mclOutput -s thaliatha obtusifolia annuum
python ../Scripts/cafe_input.py -i updated_MCLout -o genefamily/unfiltered_cafe_input.txt \
-s axillaris inflata attenuata sinuosa lycopersicum tuberosum
python ../Scripts/cafetutorial_clade_and_size_filter.py -i unfiltered_cafe_input.txt \
-o genefamily/filtered_cafe_input.txt -s

## run CAFE
python ../../../Softwares/CAFE/cafe/caferror.py -i cafe.sh -d caferror_files -v 0 -f 1
cafe caferror.sh

python python_scripts/cafetutorial_report_analysis.py -i caferror_files/cafe_final_report.cafe -o summary.cafe

