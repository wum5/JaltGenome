#!/bin/bash

#PBS -N PAML
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load python
module load biopython
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/paml/4.8 
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/hyphy-2.2.6/

## Pull out 1-to-1 orthologous sequences into gene clusters (excluding A.thaliatha and C.annumm)
cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/paml
cp ../output/pairs/mclOutput .
python ../../Scripts/remove_redundancy.py -i mclOutput -s annumm thaliatha
python ../../Scripts/1-to-1_orthologs.py -i updated_MCLout -s ../all_cds.fasta -n 6 -o gene_seqs

for i in {6001..6914};
	do prank -d=gene_seqs/gene_ortho_$i -o=from_prank/gene_ortho_$i -codon 
done

## process the sequence alignment
python ../../Scripts/rename_seqs.py -i from_prank/ -o renamed_align/
python ../../Scripts/SlidingWindows.py -i renamed_align/ -o trimed_align/ -d 6
python ../../Scripts/SlidingWindows.py -i trimed_align/ -o trimed_align2/ -d 4 -w 9
python ../../Scripts/orf_aln_process.py -i trimed_align2/ -o final_align/ -d 6

## reconstruct the gene trees
python ../../Scripts/raxml_wrapper.py --inDIR /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/paml/final_align \
--outDIR /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/paml/gene_trees \
--outgroup P.inflata,P.axillaris --cpus 2 --type DNA

## pull out the genes that have the same tree topology as the species tree
python ../../Scripts/compare_topology.py gene_trees species_tree.txt > consistent_gene_list.txt
for line in $(cat consistent_gene_list.txt); do cp final_align/${line} filter_align/${line}; done

## rename the species id 
cd filter_align
for file in *.fas; do sed -i 's/S\.lycopersicum/Slyc/g' $file; done &
for file in *.fas; do sed -i 's/S\.tuberosum/Stub/g' $file; done &
for file in *.fas; do sed -i 's/P\.inflata/Pinf/g' $file; done &
for file in *.fas; do sed -i 's/P\.axillaris/Paxi/g' $file; done &
for file in *.fas; do sed -i 's/N\.attenuata/Natt/g' $file; done &
for file in *.fas; do sed -i 's/J\.sinuosa/Jsin/g' $file; done &

## rename file name
sh rename_file.txt
cd ..

## convert to phylip format
cd ..
python ../../Scripts/seqformat_converter.py -i filter_align -o paml_input -f1 .fas -f2 .phy

## run paml (split the input data into subfolder and then use job array (#PBS -t 1-10) to submit)
cd paml_input
mkdir ch${PBS_ARRAYID}
mv Solyc${PBS_ARRAYID}g* ch${PBS_ARRAYID}/
cd ch${PBS_ARRAYID}
python ../../../../Scripts/paml_BS.py -i . -o ../paml_ch${PBS_ARRAYID}_out -t ../../species_tree.txt

## summarize the output of paml, add FDR and function annotations
for file in paml_ch*; do cat $file > paml_out
python ../../../Scripts/multiple_testing_correction.py -i paml_out -o paml_summary.txt -f function_list.txt
cd ..

## Now we want to check the effect of multinucleotide mutation (MNM) effect on selection inference
## here infileList contains the ids of all positive selected genes (pval<0.01 and FDR<0.2)
while read line; do
  cp "filter_align/"${line} "BS_MNM_test/"
done <inputList

for file in Jalsin*;
	do sed -i 's/>/#/g' $file;
	echo -e "\n\n((Paxi, Pinf), (Natt, (Jsin, (Slyc, Stub))));\n" >> $file;
done

python ../../Scripts/paml_BS_MNM.py -i /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/paml/BS_MNM_test \
-d /N/dc2/projects/jaltomt/Softwares/MNM_SelectionTests -o paml_MNM_out





