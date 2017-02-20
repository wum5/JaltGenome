# JaltGenome

## De novo asseemble genome
##### Trim low-quality Illumina short-reads
```
qsub illumina_trim.sh
```
##### Estimated genome size using k-mer
```
qsub jellyfish.sh
```
##### Assemble genome using DBG2OLC approach
```
qsub sparse_assemble.sh
qsub dbg2olc.sh
qsub blasr.sh
```
