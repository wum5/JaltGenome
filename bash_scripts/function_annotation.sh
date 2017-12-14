#!/bin/bash

#PBS -N AHRD
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=32gb
#PBS -m bea


module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
module load /N/soft/rhel6/modules/karst/DEVELOPMENT/java/jre/1.8.0_73
module load java
module unload gcc
module load gcc
module load python
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/emboss/6.5.7
module load /N/soft/rhel6/modules/mason/LIFE-SCIENCES/hmmer/3.1b2
cd /N/dc2/projects/jaltomt/GenomeAssembly/FunctionAnnotation
OD=/N/dc2/projects/jaltomt/Softwares/interproscan-5.25-64.0

ThreadNum=16

## BLAST against Arabidopsis protein database
makeblastdb -in TAIR10_pep_21101214_updated -dbtype prot -out TAIR10_pep_21101214_updated
blastp -outfmt 6 -query representative_proteins.fasta -db TAIR10_pep_20101214_updated \
-out query_vs_TAIR10.txt -num_threads $ThreadNum

## BLAST against tomato protein database
makeblastdb -in ITAG3.2_proteins.fasta -dbtype prot -out ITAG3.2_proteins.fasta
blastp -outfmt 6 -query representative_proteins.fasta -db ITAG3.2_proteins.fasta \
-out query_vs_ITAG3.2.txt -num_threads $ThreadNum

## BLAST against Sprot protein database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot.fasta
blastp -outfmt 6 -query representative_proteins.fasta -db uniprot_sprot.fasta \
-out query_vs_sprot.txt -num_threads $ThreadNum

## split input into multiple (8) files for parallel
perl ../Scripts/fasta_tools.pl --chunks 8 representative_proteins.fasta

## BLAST against Trembl protein database (each subset takes ~24h using 8 CPUs)
makeblastdb -in uniprot_trembl.fasta -dbtype prot -out uniprot_trembl.fasta
blastp -outfmt 6 -query representative_proteins_0.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_0.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_1.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_1.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_2.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_2.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_3.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_3.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_4.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_4.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_5.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_5.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_6.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_6.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_7.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_7.txt -num_threads $ThreadNum
blastp -outfmt 6 -query representative_proteins_8.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl_8.txt -num_threads $ThreadNum


## domains search using InterProscan 
## split the input files into chunks and then use job array to submit 
mkdir interproscan
perl ../Scripts/fasta_tools.pl --chunks 25 representative_proteins.fasta
cd interproscan
$OD/interproscan.sh -i representative_proteins_${PBS_ARRAYID}.fasta -f tsv -T temp${PBS_ARRAYID} 


# run AHRD to assign functional annotations
java -Xms2g -Xmx24g -jar ../../Softwares/AHRD/dist/ahrd.jar ahrd_input.yml
python ../Scripts/interpro2go.py -i ahrd_output.csv -g interpro2go > updated_ahrd_output.txt

