#!/bin/bash

#PBS -N Mapping
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load bwa
module load gcc/4.9.2 
module load star/2.5.2b
module load samtools
module load python
module load java

cd /N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/subread-1.5.3-Linux-x86_64/bin/
OD=/N/dc2/projects/jaltomt/Phylogenomics/rawdata
ThreadN=4

cd mapping
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles jalt_assembly.fa --runThreadN $ThreadN

STAR --genomeDir ./ --readFilesIn $OD/JA0702NR_shear1_p1.fastq $OD/JA0702NR_shear1_p2.fastq \
--outFileNamePrefix JA0702NR. --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
--outSAMtype SAM --runThreadN $ThreadN

STAR --genomeDir ./ --readFilesIn $OD/JA0702RP_shear1_p1.fastq $OD/JA0702RP_shear1_p2.fastq \
--outFileNamePrefix JA0702RP. --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
--outSAMtype SAM --runThreadN $ThreadN

# Converted to Bam
samtools view -bS JA0702NR.Aligned.out.sam > JA0702NR.bam
samtools view -bS JA0702RP.Aligned.out.sam > JA0702RP.bam

# Samtools sort
samtools sort -n -m 5000000000 -@ $ThreadN JA0702NR.bam JA0702NR.sort
samtools sort -n -m 5000000000 -@ $ThreadN JA0702RP.bam JA0702RP.sort

# Merge the NR and RP .bam files
samtools merge -n JA0702.sort.bam JA0702NR.sort.bam JA0702RP.sort.bam

# Get mapping statistics 
samtools flagstat JA0702.sort.bam > mapping_stats.txt
featureCounts -a jalt_assembly.aed.0.45.gff3 -o gene_expression.txt -T $ThreadN JA0702.sort.bam
