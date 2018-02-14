#!/bin/bash

#PBS -N MAKER
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load maker
module load mpich2/1.4.1-gnu
module load blat
module load gcc/4.9.2
module load python
module load biopython 
export TMPDIR=/N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation/tmp
cd /N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation


######### A function to edit CONFIG file ######### 
function set_config(){
CONFIGFILE=$1
TARGET_KEY=$2
REPLACEMENT_VALUE=$3
sed -i "s/^\($TARGET_KEY\s*=\s*\).*\$/\1$REPLACEMENT_VALUE/" $CONFIGFILE;
}

######### ab initio gene prediction ##########
# I used 16 scaffolds > 2000000 bps (totalling ~40 Mbps), my file was called scaffold_min_2000000.fa
python ../Scripts/genome_stat.py -i jalt_assembly.fa -gs 1500000000 -m 2000000

######### Edit the file so the following applies in maker_opts.ctl:
#genome=scaffold_min_2000000.fa
set_config maker_opts.ctl "genome" "scaffold_min_2000000.fa"
#est=JA0702_trinity.fasta
set_config maker_opts.ctl "est" "JA0702_trinity.fasta"
#protein=uniprot_solanaceae.fasta
set_config maker_opts.ctl "protein" "uniprot_solanaceae.fasta"
#rmlib=allRepeats.lib
set_config maker_opts.ctl "rmlib" "allRepeats.lib"
#repeat_protein=te_proteins.fasta
set_config maker_opts.ctl "repeat_protein" "te_proteins.fasta"
#snaphmm=jaltomata.cegmasnap.hmm
set_config maker_opts.ctl "snaphmm" "jaltomata.cegmasnap.hmm"
#gmhmm=gmhmm.mod
set_config maker_opts.ctl "gmhmm" "gmhmm.mod"
#est2genome=1
set_config maker_opts.ctl "est2genome" "1"
#protein2genome=1
set_config maker_opts.ctl "protein2genome" "1"
#keep_preds=1
set_config maker_opts.ctl "keep_preds" "1"
#single_exon=1
set_config maker_opts.ctl "single_exon" "1"

########### 1st round training (takes ~5h; 16 CPU) ############
mpiexec -n 16 maker --ignore_nfs_tmp </dev/null 

cd scaffold_min_2000000.maker.output
maker2zff -n -d scaffold_min_2000000_master_datastore_index.log
fathom genome.ann genome.dna -gene-stats
fathom genome.ann genome.dna -validate
cd ..
mkdir snap
cp scaffold_min_2000000.maker.output/genome.ann scaffold_min_2000000.maker.output/genome.dna snap/
cd snap
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Pult . > Pult.hmm
cd ..


########### 2nd round prediction ############

######### Edit the file to change the following applied in maker_opts.ctl:
#snaphmm=snap/Pult.hmm
set_config maker_opts.ctl "snaphmm" "snap/Pult.hmm"
#est2genome=0
set_config maker_opts.ctl "est2genome" "0"
#protein2genome=0
set_config maker_opts.ctl "protein2genome" "0"

########### 2nd round training (take ~1h; 16 CPU) ############ 
mpiexec -n 16 maker --ignore_nfs_tmp </dev/null
cd scaffold_min_2000000.maker.output
rm genome.ann
rm genome.dna
maker2zff -n -d scaffold_min_2000000_master_datastore_index.log
cd ..
mkdir snap2
cp scaffold_min_2000000.maker.output/genome.ann scaffold_min_2000000.maker.output/genome.dna snap2/
cd snap2
fathom -categorize 1000 genome.ann genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl Pult . > Pult2.hmm


########### Convert MAKER2 GFF predictions into Augustus HMM (32h) ############
cdna=/N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation/JA0702_trinity.fasta
genome=/N/dc2/projects/jaltomt/GenomeAssembly/GeneAnnotation/scaffold_min_2000000.fa
species=jaltomata
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/augustus-3.2.3/scripts
# Pay attention here (test whether could add species Dir)
export AUGUSTUS_CONFIG_PATH=/N/dc2/projects/jaltomt/Softwares/augustus-3.2.3/config

zff2gff3.pl genome.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > genome.gff3
autoAug.pl --genome=$genome --species=$species --cdna=$cdna --trainingset=genome.gff3 -v --singleCPU --useexisting
cd ..


########### MAKER run ########### 
# starting multiple maker processes in the same directory  to get it to complete
# sooner (each process will detect each other and avoid overlapping on their work, 
# so it will still get divided up efficiently) 

######### Edit the file to change the following applied in maker_exe.ctl
#tRNAscan-SE=/N/dc2/projects/jaltomt/Softwares/tRNAscan-SE-1.3.1/tRNAscan-SE #location of trnascan executable
set_config maker_exe.ctl "tRNAscan-SE" "/N/dc2/projects/jaltomt/Softwares/tRNAscan-SE-1.3.1/tRNAscan-SE"
#snoscan=/N/dc2/projects/jaltomt/Softwares/snoscan-0.9.1/snoscan #location of snoscan executable
set_config maker_exe.ctl "snoscan" "/N/dc2/projects/jaltomt/Softwares/snoscan-0.9.1/snoscan"

######### Edit the file to change the following applied in maker_opts.ctl:
#genome=jalt_assembly.fa
set_config maker_exe.ctl "genome" "jalt_assembly.fa"
#augustus_species=jaltomata
set_config maker_exe.ctl "augustus_species" "jaltomata"
#snaphmm=snap2/Pult2.hmm
set_config maker_exe.ctl "snaphmm" "snap2/Pult2.hmm"
#trna=1
set_config maker_exe.ctl "trna" "1"
#pred_stats=1
set_config maker_exe.ctl "pred_stats" "1"
#min_protein=30
set_config maker_exe.ctl "min_protein" "30"
#alt_splice=1
set_config maker_exe.ctl "alt_splice" "1"

########### run MAKER (split into 6 parts to run using 16 CPU and each takes ~100h) ############ 
fasta_tool --chunks 6 jalt_assembly.fasta 
mpiexec -n 16 maker -g jalt_assembly_0.fasta -base jalt_assembly --ignore_nfs_tmp 
mpiexec -n 16 maker -g jalt_assembly_1.fasta -base jalt_assembly --ignore_nfs_tmp
mpiexec -n 16 maker -g jalt_assembly_2.fasta -base jalt_assembly --ignore_nfs_tmp 
mpiexec -n 16 maker -g jalt_assembly_3.fasta -base jalt_assembly --ignore_nfs_tmp 
mpiexec -n 16 maker -g jalt_assembly_4.fasta -base jalt_assembly --ignore_nfs_tmp 
mpiexec -n 16 maker -g jalt_assembly_5.fasta -base jalt_assembly --ignore_nfs_tmp 
mpiexec -n 4 maker -g jalt_assembly_6.fasta -base jalt_assembly -ignore_nfs_tmp
maker -g jalt_assembly_5.fasta -base jalt_assembly --ignore_nfs_tmp -tries 5
maker -dsindex -g jalt_assembly.fasta --ignore_nfs_tmp

cd jalt_assembly.maker.output
gff3_merge -d jalt_assembly_master_datastore_index.log
fasta_merge -d jalt_assembly_master_datastore_index.log

########### extract AED < 1 genes only from MAKER output ########### 
PATH=$PATH:/N/dc2/projects/jaltomt/GenomeAssembly/Scripts
## extract genes with AED cutoff of 0.6 to get a comparable number of genes relative other Solanaceae genome
python ../../Scripts/annotation_by_aed.py --gff jalt_assembly.all.gff --AED 0.6 --out
## get a summary of gene annotation
python ../../Scripts/annotation_stats.py --gff jalt_assembly.aed.0.6.gff3 

## extract the representative (longest) protein sequences 
perl -lne 'print $1 if /\tmRNA\t.+ID=([^;]+).+_AED=(.+?);/ and $2<=0.6' jalt_assembly.all.gff > jalt_assembly.aed-0.6.ids
perl ../../Scripts/fastaqual_select.pl -f jalt_assembly.all.maker.proteins.fasta -inc jalt_assembly.aed-0.6.ids > jalt_assembly.proteins.aed-0.6.fasta
python ../../Scripts/representative_transcripts.py --fasta jalt_assembly.proteins.aed-0.6.fasta

## change the gene IDs into more human readable names 
python ../../Scripts/convert_geneids.py -i representative_proteins.fa -g jalt_assembly.aed.0.6.gff3

## extract cds for representative genes
../../Softwares/gffread/gffread -x all_cds -g jalt_assembly.fa updated_geneids.gff
python ../Scripts/representative_transcripts.py --fasta all_cds.fa --type cds


