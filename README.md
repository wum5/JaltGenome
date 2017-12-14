# JaltGenome

## Overview
* Raw scripts/Pipeline for the "Jaltomata Genome" Project.
* detailed information for each step are recorded in the corresponding bash script
* Each bash script might combine scripts for multiple runs (need to check before running)
* Still in updating!

## Contributors 
* Meng Wu
* https://github.com/wum5/JaltGenome

## De novo asseemble genome
##### Assemble genome using Masurca approach
```
qsub masurca.sh
```
##### Genome assembly evaluation (conserved single-copy orthologs and RNA-seq uniquely mapped)
```
python genome_stat.py -i final_assembly.fasta -s 1500000000
qsub cegs_eval.sh
qsub mapping.sh
```
##### Remove contaminants in the assembly and assemble Organelle genome separately
```
qsub blast_contaminants.sh
qsub organelle_assembly.sh
```
##### RepeatMasking and train gene predictors
```
qsub repeat_annot2.sh
```
##### Genome annotation using the pipeline Maker (need to change maker_opts.ctl step by step; detail in bash script)
```
qsub maker.sh
```
##### Functional Annotations
```
qsub function_annotation.sh
```
##### Phylogenetic analyses
```
qsub ortholog_inference.sh
qsub phylogeny_reconstruction.sh
```
##### Gene family analyses
```
qsub gene_family_analyses.sh
qsub go_enrichment.sh
##### Examination on SEUSS gene
```
qsub candidate_analyze.sh
```
#####  Molecular evolution (dN/dS) analyses
```
qsub dnds_test.sh
```

