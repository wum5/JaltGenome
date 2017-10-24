#!/bin/bash

#PBS -N REPEAT
#PBS -l nodes=1:ppn=16,walltime=96:00:00,vmem=64gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load gcc/4.9.2
module load perl
module load /N/soft/rhel6/modules/mason/LIFE-SCIENCES/bioperl
module load /N/soft/rhel6/modules/mason/LIFE-SCIENCES/maker/2.31.6
module load /N/soft/rhel6/modules/mason/LIFE-SCIENCES/hmmer/3.0
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/muscle/3.8.31
module load /N/soft/rhel6/modules/karst/LIFE-SCIENCES/ncbiblast+/2.2.31

seqfile='jalt_assembly.fa'
seqfileindex='jalt_assembly'


### This pipeline was modified based on the advanced repeat masking steps provided by MAKER-P
### MITE Hunter: take ~120h using 16 CPU 
cd /N/dc2/projects/jaltomt/GenomeAssembly/RepeatAnnotation/MITE
OD=/N/dc2/projects/jaltomt/Softwares/MITE_Hunter

#perl $OD/MITE_Hunter_manager.pl -i ../$seqfile -g jaltomata -n 16 -S 12345678 
DIR_CRL=/N/dc2/projects/jaltomt/GenomeAssembly/RepeatAnnotation/CRL_Scripts1.0
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/RepeatModeler


## LTR Harvest (young LTRs): take ~30h 
## This step is replaced by LTR_retriver now...

### Collecting repetitive sequences by RepeatModeler: take ~300h 
cd ..
#mkdir updated_mask
cd updated_mask
cp ../MITE/Mite.lib .
cp ../jalt_assembly.fa .
cp ../classification/jaltomata/jalt_assembly.fa.mod.masked .
cp ../classification/jaltomata/jalt_assembly.fa.mod.LTRlib.fa ./allLTR.lib

RepeatMasker -lib MITE.lib -dir . jalt_assembly.fa.mod.masked -pa 15

perl $DIR_CRL/rmaskedpart.pl $seqfile.masked 50  >  umseqfile

## Use 16 CPU to speed up RepeatModeler processing
BuildDatabase -name umseqfiledb -engine ncbi umseqfile
RepeatModeler -database umseqfiledb -pa 15 -engine ncbi >& umseqfile.out
RepeatClassifier -consensi consensi.fa  ## this step taks a long time and the input can be splited into multiple files


cp RM_43234.SatSep232036122017/consensi.fa.classified ./
perl $DIR_CRL/repeatmodeler_parse.pl --fastafile consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta  \
--identities repeatmodeler_identities.fasta 

makeblastdb -in Tpases020812DNA -dbtype prot
blastx -query repeatmodeler_unknowns.fasta -db Tpases020812DNA -evalue 1e-10 -num_descriptions 10 \
-out modelerunknown_blast_results.txt

perl $DIR_CRL/transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown \
repeatmodeler_unknowns.fasta

mv unknown_elements.txt ModelerUnknown.lib
cat identified_elements.txt repeatmodeler_identities.fasta > ModelerID.lib
cat MITE.lib allLTR.lib > allMITE_LTR.lib
cat allMITE_LTR.lib ModelerID.lib > KnownRepeats.lib

### Exclusion of gene fragments
makeblastdb -in alluniRefprexp070416 -dbtype prot
blastx -query ModelerUnknown.lib -db alluniRefprexp070416 -evalue 1e-10 -num_descriptions 10 \
-out ModelerUnknown.lib_blast_results.txt
blastx -query KnownRepeats.lib -db alluniRefprexp070416 -evalue 1e-10 -num_descriptions 10 \
-out KnownRepeats.lib_blast_results.txt

ProtExcluder1.2/ProtExcluder.pl -option ModelerUnknown.lib_blast_results.txt ModelerUnknown.lib
ProtExcluder1.2/ProtExcluder.pl -option KnownRepeats.lib_blast_results.txt KnownRepeats.lib

cat KnownRepeats.lib ModelerUnknown.lib > allRepeats.lib


### RepeatMask all repeats: take ~300h 
RepeatMasker -lib allRepeats.lib -dir . -pa 16 -gff -excln jalt_assembly.fa 

PATH=$PATH:/N/dc2/projects/jaltomt/GenomeAssembly/RepeatAnnotation
faToTwoBit jalt_assembly.fa.mod jalt_assembly.2bit
perl /N/dc2/projects/jaltomt/Softwares/RepeatMasker/util/buildSummary.pl -species jaltomata \
-genome jalt_assembly.2bit jalt_assembly.fa.out > jalt_assembly.fa.tbl


