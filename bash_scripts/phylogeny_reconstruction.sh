#!/bin/bash

#PBS -N Phylogeny
#PBS -l nodes=1:ppn=1,walltime=48:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

set -e
set -u
set -o pipefail

module load python
module load biopython
module load raxml
module load java


cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/phylogeny
OS=/N/dc2/projects/jaltomt/Softwares/ASTRAL
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/mrbayes-3.2.6/src
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/bucky-1.4.4/src


## Concatenated the sequences
python concatenate_matrix.py renamed_align 100 8 aa ../phylogeny/concatenation
python ../../Scripts/rename_seqs.py -i trimed_align/ -o renamed_align/


## Generate the RAxML phylogeny
raxmlHPC-PTHREADS -T 8 -f a -x 12345 -N 100 -p 12345 -m GTRGAMMA -q concatenated.model \
-s concatenated.phy -o A.thaliatha -n supermatrix


## Generate gene trees with bootstrap
python ../../Scripts/raxml_wrapper.py --inDIR /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/alignments/final_align/ \
--outDIR /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/phylogeny/genetree/ --outgroup A.thaliatha \
--cpus 8 --type DNA --bootstrap


## Generate majority consensus tree
cat genetree/* > genetrees.tre
raxmlHPC -L MRE -z genetrees.tre -m GTRGAMMA -n T1

## parse the proportions of gene trees supported for the ambiguous internode
python ../../Scripts/parse_internode.py -i genetrees.tre -s C.annuumm J.sinuosa S.lycopersicum S.tuberosum


## Generate ASTRAL quartet-based trees
cd genetree
# The genetree folder contains all the best RAxML gene tree files for each gene
find ./ -name "*bestTree*" | sort -V | xargs cat > ../genetrees.tre
cd ..
# The bootstrap folder contains all the bootstraps files generated for each gene trees
find bootstrap/ -name "*bootstrap*" | sort -V > bsfile.txt

java -Xincgc -Xms3000m -Xmx3000m -XX:MaxPermSize=3000m -jar $OS/astral.4.10.12.jar  \
-i genetrees.tre -b bsfile.txt -o astral_out.tre -r 100


## Generate Bucky concordance tree
outgroup="A.thaliatha"

## create input files for MrBayes
cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs
python ../Scripts/seqformat_converter.py -i alignments/final_align/ -o phylogeny/bucky/ -f1 .fas -f2 .nex

cd phylogeny/bucky/
text="\nbegin mrbayes;\n\toutgroup ${outgroup};\n\tset autoclose=yes nowarnings=yes;\n\tlset nst=6 rates=invgamma; [specifies a GTR+I+G]\n\tmcmc ngen=1000000 samplefreq=1000 printfreq=1000 savebrlens=yes filename=mymb;\n\tquit;\nend;"
for file in *.nex; do echo -e $text >> $file; done
python directory_subpackage.py ./ 1 .nex

## run MrBayes (can be splited into subset of genes for parallel)
num=$(ls -d sub_folder_* | wc -l)
for d in sub_folder_{1..$num}; do cd $d; mb *.nex; cd ..; done

## run BUCKy 
for d in sub_folder{1..3441}; do echo $d"/mymb" >> bucky_infilelist; done
for d in sub_folder*; do cd $d; mbsum -n 501 -o mymb mymb.run?.t; cd ..; done
bucky -n 1000000 -sg -i bucky_infilelist



