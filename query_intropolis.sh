#!/bin/sh
#$ -N pandey_intropolis
#$ -l mem_free=10G,h_vmem=11G,h_fsize=50G
#$ -m e
#$ -pe local 12
#$ -cwd

## Create logs dir
mkdir -p logs

echo "**** Job starts ****"
date

## Create jx.tsv file
module load R/3.3.x
Rscript query_intropolis.R

## Query intropolis
cat <(gzip -cd /dcl01/leek/data/sra_work/v2/intropolis.v2.hg38.tsv.gz | cut -f1-4) jx.tsv | /users/lcollado/software/coreutils/bin/sort --parallel=12 -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c | awk '$1 == 2' > jx_in_intropolis.tsv

echo "**** Job ends ****"
date

## Move log files
mv pandey_intropolis.* logs/
