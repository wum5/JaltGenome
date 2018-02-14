#!/bin/bash

#PBS -N Hetero
#PBS -l nodes=1:ppn=4,walltime=12:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/samtools/0.1.18
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/vcftools/0.1.10
module unload python
module load python/3.5.0
OS=/N/dc2/projects/jaltomt/Softwares/mvftools
REF=/N/dc2/projects/jaltomt/GenomeAssembly/KEY_FILES/jalt_assembly.fa

## Using RNA-seq data
cd /N/dc2/projects/jaltomt/GenomeAssembly/Mapping/rna/HeterData

samtools mpileup -uD -f $REF JA0701.sort.bam JA0456.sort.bam JA0694.sort.bam \
JA0450.sort.bam JA0798.sort.bam JA0711.sort.bam JA0723.sort.bam JA0608.sort.bam \
JA0702.sort.bam JA0726.sort.bam JA0432.sort.bam JA0816.sort.bam JA0719.sort.bam \
JA0010.sort.bam | bcftools view -cg - > var.vcf

## split the huge file into multiple smaller files
head -n 10000 var.vcf | grep "^#" > header
grep -v "^#" var.vcf > variants
split -l 100000000 variants
for i in x*; do cat header $i >$i.vcf && rm -f $i; done
rm -f header variants

## transform from vcf to mvf
filename=xaj
python3 $OS/vcf2mvf.py --vcf ${filename}.vcf --out ${filename}.mvf --lowdepth 10 --lowqual 30

## do the calculation of Heterozygosity
python ../../../Scripts/ancestral_variation.py -i var.mvf -t species_hetero


## Using DNA-seq data
cd /N/dc2/projects/jaltomt/GenomeAssembly/Mapping/dna
samtools view -s 3.25 -@ 8 -b dna_seqs.sort.bam > 30X_sort.bam
INFILE=30X_sort.bam
samtools mpileup -uD -Q 20 -q 20 -f $REF $INFILE | bcftools view -cg -N - > 30X_var.vcf

## estimate the base call error rate
mkdir temp
cd temp
grep -v '^#' ../30X_var.vcf | split -l 10000000
for file in *; do cat ../header.txt ${file} > ${file}".vcf"; done

find . -name "*.var" | xargs basename -s ".var" | \
xargs -P 4 -I{} vcftools --vcf {}.var --min-meanDP 5 --recode --out {}_DP5

INFILE=all.recode.vcf
## total sites: 1373175316
find . -name "*_DP5.recode.vcf" | xargs -P 32 grep -v '^#' | wc -l
## heterozygous sites: 321364
find . -name "*_DP5.recode.vcf" | xargs -P 32 grep '0/1:' | wc -l
## homozygous alterative sites: 41185
find . -name "*_DP5.recode.vcf" | xargs -P 32 grep '1/1:' | wc -l
## error rate calculation: 41185.0/(1373175316-321364)*100=0.0029999


