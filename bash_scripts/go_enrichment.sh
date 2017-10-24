#!/bin/bash

#PBS -N GO_ENRICHMENT
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


OD=/N/dc2/projects/jaltomt/Softwares/ontologizer
cd /N/dc2/projects/jaltomt/GenomeAssembly/Orthologs/

# prepare go mapping file
cut -f1,6 ../KEY_FILES/ahrd_updated.txt > go_test/jalt_sin_go.ids

## pull out sinuosa-specific gene families
python ../Scripts/lineage_specific_genes.py -i genefamily/updated_MCLout -s sinuosa

## pull out genes in significantly expanded/contracted families in sinuosa
python ../../Scripts/fast_evol_family.py -i expanded_families_jalt.txt \
-c updated_MCLout -g ahrd_updated.txt -s sinuosa -o expanded_genes.txt
python ../../Scripts/fast_evol_family.py -i contracted_families_jalt.txt \
--c updated_MCLout --g ahrd_updated.txt -s sinuosa -o contracted_genes.txt

## prepare inputfile for Ontologizer
cut -f2,7 contracted_genes.txt > contracted_go.txt
cut -f2,7 expanded_genes.txt > expanded_go.txt

## Run Ontologizer v2.0
java -Xmx258m -jar $OD/Ontologizer.jar -a jalt_sin_go.ids -g go.obo -s expanded_go.txt \
-p pop_genes.txt -c Parent-Child-Union -m Benjamini-Hochberg
java -Xmx258m -jar $OD/Ontologizer.jar -a jalt_sin_go.ids -g go.obo -s contracted_go.txt \
-p pop_genes.txt -c Parent-Child-Union -m Benjamini-Hochberg
java -Xmx258m -jar $OD/Ontologizer.jar -a jalt_sin_go.ids -g go.obo -s sinuosa_specific.txt \
-p pop_genes.txt -c Parent-Child-Union -m Benjamini-Hochberg