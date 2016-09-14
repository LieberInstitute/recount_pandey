#!/bin/sh

## Usage:
# bash recount_pandey.sh

## Example in JHPCE
module load R/3.3.x
module load ucsctools
module load wiggletools/default

## Download required scripts
wget https://raw.githubusercontent.com/leekgroup/recount-website/5881e1359fb470c529a011e5440323d5247f7ca6/recount-prep/prep_merge.R .
wget https://raw.githubusercontent.com/leekgroup/recount-website/5881e1359fb470c529a011e5440323d5247f7ca6/recount-prep/prep_sample.R .
wget https://raw.githubusercontent.com/leekgroup/recount-website/5881e1359fb470c529a011e5440323d5247f7ca6/recount-prep/prep_setup.R .

## Some common variables
DATADIR="/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs"
COUNTS="/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/counts.tsv.gz"
BWTOOL="/dcl01/leek/data/bwtool/bwtool-1.0/bwtool"
WIGGLE="wiggletools"

## Download some required files
Rscript prep_setup.R

## Display help info on how to run prep_sample.R
Rscript prep_sample.R -h

## Process all samples
cut -f 5 /dcl01/leek/data/sunghee/all_s3.manifest | parallel --jobs 5 Rscript prep_sample.R -f ${DATADIR}/{}.bw -c ${COUNTS} -b ${BWTOOL} -w ${WIGGLE} -a TRUE

## Now merge results
paste rse_temp/counts_exon_* > counts_exon.tsv
gzip counts_exon.tsv
paste rse_temp/counts_gene_* > counts_gene.tsv
gzip counts_gene.tsv

## Display help info on how to run prep_merge.R
Rscript prep_merge.R -h

## Merge rse objects and create junction rse object
BWDIR="/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs"
JUNCTIONS="/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/first_pass_junctions.tsv.gz"
MANIFEST="/dcl01/leek/data/sunghee/all_s3.manifest"
WIGTOBIGWIG="wigToBigWig"
Rscript prep_merge.R -b ${BWDIR} -j ${JUNCTIONS} -m ${MANIFEST} -w ${WIGGLE} -t ${WIGTOBIGWIG} -c TRUE

## Clean up external scripts
rm prep_merge.R prep_sample.R prep_setup.R
