#!/bin/bash

#PBS -N CandiGenes
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31
module load samtools
module load raxml
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/prank/bin

cd /N/dc2/projects/jaltomt/GenomeAssembly/candidate_genes

filename="suess_all"
## prepare sequence alignment
prank -d=${filename}".fas" -o=${filename}".aln" -DNA 

cd ../../Softwares/phyutility
OD=/N/dc2/projects/jaltomt/GenomeAssembly/candidate_genes
./phyutility -clean 0.25 -in $OD/${filename}".aln.best.fas" -out $OD/${filename}".aln-pht"
cd /N/dc2/projects/jaltomt/GenomeAssembly/candidate_genes
raxmlHPC -m GTRGAMMA -n ${filename} -s ${filename}".aln-pht" -p 12345  -o Peaxi162Scf00183g00037 Peinf101Scf00652g10019


# Blast to PacBio reads 
makeblastdb -in ../../PacBioReads/PacBio.fasta -dbtype nucl -out PacBio.fasta
blastn -outfmt 6 -query SEU_tomato.fa -db PacBio.fasta -out scf29960_blast.txt -task dc-megablast


## check the read depth of investigated regions to see whether there are depth peaks 
## corresponding to the duplicated regions
samtools depth dna_seqs.sort.bam -r scf7180000031961 > scf31961_read_depth.txt
samtools depth dna_seqs.sort.bam -r scf7180000029960 > scf29960_read_depth.txt

python ../Scripts/read_depth_plot.py -i scf31961_read_depth.txt -o scf31961_depth.jpg
python ../Scripts/read_depth_plot.py -i scf29960_read_depth.txt -o scf29960_depth.jpg

## run Rscript "seuss_plot.R"

