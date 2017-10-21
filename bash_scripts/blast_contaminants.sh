#!/bin/bash

#PBS -N BLAST-CONTAM
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
cd /N/dc2/projects/jaltomt/GenomeAssembly/BLAST_ALL

makeblastdb -in nt -parse_seqids -dbtype nucl

awk 'BEGIN {n_seq=0;} /^>/ {if(n_seq%500==0){file=sprintf("myseq%d.fa",n_seq);} \
print >> file; n_seq++; next;} { print >> file; }' < init_assembly.fa

blastn -task megablast -num_threads 8 -query myseq0.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out1
blastn -task megablast -num_threads 8 -query myseq6000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out13

blastn -task megablast -num_threads 8 -query myseq500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out2
blastn -task megablast -num_threads 8 -query myseq5500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out12

blastn -task megablast -num_threads 8 -query myseq6500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out14
blastn -task megablast -num_threads 8 -query myseq7000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out15

blastn -task megablast -num_threads 8 -query myseq1000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out3
blastn -task megablast -num_threads 8 -query myseq5000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out11

blastn -task megablast -num_threads 8 -query myseq1500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out4
blastn -task megablast -num_threads 8 -query myseq4500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out10

blastn -task megablast -num_threads 8 -query myseq2000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out5
blastn -task megablast -num_threads 8 -query myseq4000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out9

blastn -task megablast -num_threads 8 -query myseq2500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out6
blastn -task megablast -num_threads 8 -query myseq3500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out8

blastn -task megablast -num_threads 8 -query myseq3000.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out7
blastn -task megablast -num_threads 8 -query myseq7500.fa -db nt -outfmt 6 -evalue 1e-10 -out blastn.out16

### extract the best blast hit based on combined hit score with the gi number
cat blastn.out* > blast_out.txt
python ../Scripts/combine_blastout.py blast_out.txt

## Search the rest ones (not found by megablast) by blastn 
python ../Scripts/remove_seqs.py init_assembly.fa remove_seqs.txt rest.fa
blastn -task blastn -query rest.fa -db nt -outfmt 6 -evalue 1e-5 -out blastn.out18
python ../Scripts/combine_blastout.py blastn.out17


## Remove contaminants
python ../Scripts/remove_seqs.py init_assembly.fa remove_seqs.txt jalt_assembly_rm.fa

