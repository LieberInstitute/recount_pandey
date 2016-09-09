#!/bin/sh
#$ -N pandey_recount
#$ -l mem_free=250G,h_vmem=300G,h_fsize=40G
#$ -m e
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
