#!/bin/bash

#PBS -N Mapping
#PBS -l nodes=1:ppn=8,walltime=8:00:00,vmem=32gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load bwa
module load star
module load samtools
module load python
module load java
module load trimmomatic


PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/subread-1.5.3-Linux-x86_64/bin/
#OD=/N/dc2/projects/jaltomt/Phylogenomics/rawdata
OD=/N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation
ThreadN=8


### DNA mapping 
#cd /N/dc2/projects/jaltomt/GenomeAssembly/IlluminaReads
trimmomatic-0.36 PE -phred33 -threads $ThreadN raw_reads1.fq raw_reads2.fq xfiltered_reads1_paired.fq \
xfiltered_reads1_unpaired.fq xfiltered_reads2_paired.fq xfiltered_reads2_unpaired.fq \
ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

cd ../Mapping/dna 
OD=/N/dc2/projects/jaltomt/GenomeAssembly/IlluminaReads
bwa index jalt_assembly.fa
bwa mem jalt_assembly.fa $OD/xfiltered_reads1_paired.fq $OD/xfiltered_reads2_paired.fq -t $ThreadN > dna_seqs.sam
samtools view -bS -@ $ThreadN dna_seqs.sam | samtools sort -m 5000000000 > dna_seqs.sort.bam
samtools flagstat dna_seqs.sort.bam > mapping_stats.txt
samtools index dna_seqs.sort.bam

## extract uniquely mapped paired reads
samtools view -h dna_seqs.sort.bam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -q 10 \
-F 3332 -f 2 -b > dna_seqs_sorted_properlypaired.bam


### RNA	mapping
cd /N/dc2/projects/jaltomt/GenomeAssembly/Mapping/rna
OD=/N/dc2/projects/jaltomt/Phylogenomics/KEY_FILES/rawdata

#STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles jalt_assembly.fa --runThreadN $ThreadN

## Mapping different libraries including:JA0432,JA0450,JA0456,JA0608,JA0694,JA0701,JA0702,
## JA0711,JA0719,JA0723,JA0726,JA0798,JA0816
sample_names=(JA0432 JA0450 JA0456 JA0608 JA0694 JA0701 JA0702 JA0711 JA0719 JA0723 JA0726 JA0798 JA0816)

for index in ${sample_names[@]}
do
	STAR --genomeDir ./ --readFilesIn $OD/$index'NR_shear1_p1.fastq' $OD/$index'NR_shear1_p2.fastq' \
	--outFileNamePrefix $index'NR'. --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
	--outSAMtype SAM --runThreadN $ThreadN

	STAR --genomeDir ./ --readFilesIn $OD/$index'RP_shear1_p1.fastq' $OD/$index'RP_shear1_p2.fastq' \
	--outFileNamePrefix $index'RP.' --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
	--outSAMtype SAM --runThreadN $ThreadN

	# Converted to Bam and then sort
	samtools view -bS -@ $ThreadN $index'NR.Aligned.out.sam' | samtools sort -m 5000000000 > $index'NR.sort.bam'
	samtools view -bS -@ $ThreadN $index'RP.Aligned.out.sam' | samtools sort -m 5000000000 > $index'RP.sort.bam'

	# Merge the NR and RP .bam files
	samtools merge -@ $ThreadN $index'.sort.bam' $index'NR.sort.bam' $index'RP.sort.bam'
	mv $index'.sort.bam' finalData/
done


# Get mapping statistics 
cd finalData
featureCounts -a ../updated_geneids.gff3 -o uniquely_mapped_counts.txt \
-g ID -t gene -T $ThreadN *.sort.bam
featureCounts -a ../updated_geneids.gff3 -o multiply_mapped_counts.txt \
-M -g ID -t gene -T $ThreadN *.sort.bam

