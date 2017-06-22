# JaltGenome

## De novo asseemble genome
##### Assemble genome using DBG2OLC approach
```
sh illumina_trim.sh
sh sparse_assemble.sh
sh dbg2olc.sh
sh blasr.sh
```
##### Assemble genome using Masurca approach
```
sh masurca.sh
```
##### Genome assembly evaluation
```
python genome_stat.py -i final_assembly.fasta -s 1500000000
sh cegs_eval.sh
sh mapping.sh
```
##### Remove contaminants in the assembly and assemble Organelle genome separately
```
sh blast_contaminants.sh
sh organelle_assembly.sh
```
##### Genome annotation using the pipeline Maker (need to change maker_opts.ctl step by step; detail in bash script)
```
sh repeat_annot.sh
sh ab_initio_traning.sh
sh maker.sh
```
