Recount-style output for Pandey's data
======================================

This document describes the files at /dcl01/leek/data/recount_pandey created by Leonardo Collado Torres. Please get in touch with him if you have any questions regarding these files.

## Recount-style output

The scripts `main_job.sh` and `recount_pandey.sh` create objects similar to the ones from the [`recount` project](https://jhubiostatistics.shinyapps.io/recount/) using scripts available at [recount-website/recount-prep](https://github.com/leekgroup/recount-website/tree/master/recount-prep). The files that get created are:

```
counts_exon.tsv.gz
counts_gene.tsv.gz
rse_exon.Rdata
rse_gene.Rdata
rse_jx.Rdata
bw
    mean.bw
    
## Log files
logs/pandey_recount*
```

## Which exon-exon junctions overlap Intropolis v2

The scripts `query_intropolis.sh`, `query_intropolis.R` and `classify_intropolis.R` are all used together in order to check which exon-exon junctions from the `rse_jx.Rdata` file are present in _Intropolis v2_. The main output of this script is `rse_jx_with_intropolis.Rdata` which is just like `rse_jx.Rdata` but includes the logical column `is_intropolis` labeling if the exon-exon junction was found in Intropolis v2.

The full list of outputs from these scripts are:

```
## Temporary files used by the scripts
jx_in_intropolis.tsv.gz
jx.tsv.gz

## Log files
logs/prandey_intropolis*

## Final output
rse_jx_with_intropolis.Rdata
```

## Exploring junction results

The R script `characterize_intropolis.R` explores the file `rse_jx_with_intropolis.Rdata` and creates a series of exploratory plots saved in three different PDF files. These are:

```
exploratory_plots.pdf
maximum_coverage_jx_Intropolis_not_annotated_UCSC.pdf
maximum_coverage_jx_new.pdf

## Log file:
logs/characterize_intropolis_log.txt
```


## Merging by cell

The R script `merge_cells.R` uses the information from `CellMap_codes.csv` to merge the technical replicates. It creates the `rse_with_cell` directory which contains the _RangedSummarizedExperiment_ objects with the cell information but prior to merging by cell type. These files could be useful for some quality control or other analyses. Then `merge_cells.R` creates the output directory `rse_merged` with the RSE objects that have merged the information from the cells (that is, 34 columns instead of 258).

The files in `rse_merged` are:

```
rse_merged
    rse_exon.Rdata
    rse_gene.Rdata
    rse_jx.Rdata
    
## Log file:
logs/merge_cells_log.txt
```

## Exploring junction results by cells

Then, the script `rse_merged/characterize_intropolis_merged.R` used the information from `rse_merged/rse_jx.Rdata` to create another set of exploratory plots (analogous to the ones created previously) to explore the exon-exon junctions found against the number of cells (instead of _samples_ or technical replicates). It creates the PDF files:

```
rse_merge
    exploratory_plots_merged.pdf
    maximum_coverage_jx_Intropolis_not_annotated_UCSC_merged.pdf
    maximum_coverage_jx_new_merged.pdf

## Log file:
rse_merged/logs/characterize_intropolis_merged_log.txt
```

