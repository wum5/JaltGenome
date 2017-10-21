#!/bin/bash

#PBS -N AHRD
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=160gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
module load /N/soft/rhel6/modules/karst/DEVELOPMENT/java/jre/1.8.0_73
module load java
module unload gcc
module load gcc/4.9.2
module load python
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/emboss/6.5.7
module load /N/soft/rhel6/modules/mason/LIFE-SCIENCES/hmmer/3.1b2
cd /N/dc2/projects/jaltomt/GenomeAssembly/FunctionAnnotation
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/interproscan-5.25-64.0

ThreadNum=8

# BLAST against Arabidopsis protein database
makeblastdb -in TAIR10_pep_21101214_updated -dbtype prot -out TAIR10_pep_21101214_updated
blastp -outfmt 6 -query representative_proteins.fasta -db TAIR10_pep_21101214_updated \
-out query_vs_TAIR10.txt -num_threads $ThreadNum

# BLAST against tomato protein database
makeblastdb -in ITAG2.4_proteins.fasta -dbtype prot -out ITAG2.4_proteins.fasta
blastp -outfmt 6 -query representative_proteins.fasta -db ITAG2.4_proteins.fasta \
-out query_vs_ITAG2.4.txt -num_threads $ThreadNum

# BLAST against Sprot protein database
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot.fasta
blastp -outfmt 6 -query representative_proteins.fasta -db uniprot_sprot.fasta \
-out query_vs_sprot.txt -num_threads $ThreadNum

# BLAST against Trembl protein database
makeblastdb -in uniprot_trembl.fasta -dbtype prot -out uniprot_trembl.fasta
blastp -outfmt 6 -query representative_proteins.fasta -db uniprot_trembl.fasta \
-out query_vs_trembl.txt -num_threads $ThreadNum

# domains search using InterProscan
interproscan.sh -i representative_proteins.fasta -f tsv

# run AHRD to assign functional annotations
java -Xms24g -Xmx120g -jar ../../Softwares/AHRD/dist/ahrd.jar ahrd_input.yml
