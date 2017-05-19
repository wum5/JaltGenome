# JaltGenome

## De novo asseemble genome
##### Assemble genome using DBG2OLC approach
```
qsub illumina_trim.sh
qsub sparse_assemble.sh
qsub dbg2olc.sh
qsub blasr.sh
```
##### Genome assembly evaluation
```
python genome_stat.py -i final_assembly.fasta -s 1500000000
qsub cegs_eval.sh
```
##### Genome annotation using the pipeline Maker (need to change maker_opts.ctl step by step; detail in bash script)
```
qsub repeatmodeler.sh
qsub maker.sh
```
