# JaltGenome

## De novo asseemble genome
##### Trim low-quality Illumina short-reads
```
qsub illumina_trim.sh
```
##### Estimated genome size using k-mer
```
qsub jellyfish.sh
rstrip Kmer_GS.R
```
##### Assemble genome using DBG2OLC approach
```
qsub sparse_assemble.sh
qsub dbg2olc.sh
qsub blasr.sh
```
##### Genome assembly evaluation
```
python genome_stat.py -i ../../GenomeAssembly/Assembly/final_assembly.fasta -s 1400000000
qsub cegs_eval.sh
qsub bwa_map.sh
```
##### Genome annotation using the pipeline Maker (detail in bash script)
```
qsub maker.sh
```
