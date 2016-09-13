#!/bin/sh
#$ -N pandey_recount
#$ -l mem_free=10G,h_vmem=11G,h_fsize=40G
#$ -m e
#$ -pe local 10
#$ -cwd

## Create logs dir
mkdir -p logs

echo "**** Job starts ****"
date

bash recount_pandey.sh

echo "**** Job ends ****"
date

## Move log files
mv pandey_recount.* logs/
