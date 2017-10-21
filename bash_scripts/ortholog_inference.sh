#!/bin/bash

#PBS -N OrthoMCL
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load python
module load biopython

PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/mcl-14-137/bin
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/interproscan-5.25-64.0/bin/blast/2.2.24/bin

export PERL5LIB="/N/u/wum5/Mason/perl5/lib/perl5:$PERL5LIB" 
export PATH="/N/u/wum5/Mason/perl5/bin:$PATH"
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/prank/bin
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/gffread
  
  
## Installation
OD=/N/dc2/projects/jaltomt/Softwares/orthomcl-pipeline
perl $OD/scripts/orthomcl-pipeline-setup.pl
perl $OD/scripts/orthomcl-setup-database.pl --user root --password rip42ua --host rdc04.uits.iu.edu:3136 \
--database orthomcl --outfile $OD/orthomcl.conf 


## Runs orthomcl using the input fasta files under input/ and orthomcl.confg as config file.
## Places data in output/.  Gets other parameters (blast, etc) from default config file.
cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs
$OD/bin/orthomcl-pipeline -i ./input -o ./output -m $OD/orthomcl.conf --nocompliant --yes


## Extract cds for jaltomata and download other cds from solgenomics for following phylogeny analyses
gffread -w jalsin_cds.fa -g ../GeneAnnotation/jalt_assembly.fa ../GeneAnnotation/jalt_assembly.aed.0.45.gff3
python ../Scripts/representative_transcript.py --fasta jalsin_cds.fa --gff jalt_assembly.aed.0.45.gff3 \
--outfile representative_cds.fasta



## Pull out 1-to-1 orthologous sequences into gene clusters
mkdir 1-to-1
cat input/*fasta > all_input.fasta
python ../Scripts/1-to-1_orthologs.py -i output/pairs/mclOutput2 -s all_input.fasta -n 8 -o 1-to-1

## Make sequence alignment using Prank v.140603
mkdir alignment 
make alignment/from_prank
for i in {1..500};
	do prank -d=1-to-1/gene_ortho_$i -o=alignments/from_prank/gene_ortho_$i -codon 
done

## Trim alignment using 15-bp sliding window
cd alignments
python ../../Scripts/rename_seqs.py -i from_prank/ -o renamed_align/
python ../../Scripts/SlidingWindows.py -i renamed_align/ -o trimed_align/ -e A.thaliatha -d 6
python ../../Scripts/orf_aln_process.py -i trimed_align/ -o final_align/ -d 8
python concatenate_matrix.py final_align/ 200 8 dna ../phylogeny/concatenated


## Pick sequence alignments >200 amino acids after removing gaps or missing sites
mkdir trimed_align
python phyutility_wrapper.py from_prank/ trimed_align/ 0.95 aa
