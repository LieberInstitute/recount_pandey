## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading rse_jx.Rdata'))
system.time( load('rse_jx.Rdata') )

## Get junctions in the data
jx <- rowRanges(rse_jx)
jx_df <- data.frame(chr = seqnames(jx), start = start(jx), end = end(jx), strand = strand(jx))

## Write to disk
write.table(jx_df, file = 'jx.tsv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

## Reproducibility info
proc.time()
options(width = 120)
session_info()
