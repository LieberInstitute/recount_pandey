## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading rse_jx_with_intropolis.Rdata'))
system.time( load('rse_jx_with_intropolis.Rdata') )

## Keep only those not in intropolis
rse_new <- subset(rse_jx, subset = !in_intropolis)
rse_present <- subset(rse_jx, subset = in_intropolis)



n_samples <- rowSums(assays(rse_new, 1)$counts > 0)
table(n_samples)
summary(n_samples)

n_samples_p <- rowSums(assays(rse_present, 1)$counts > 0)
table(n_samples_p)
summary(n_samples_p)

## What is going on with some that have 0?
weird <- which(n_samples == 0)
weird_gr <- rowRanges(rse_new)[weird]

opt <- list(jx_file = '/dcl01/leek/data/sunghee_analysis/processed/cross_sample_results/first_pass_junctions.tsv.gz')
message(paste(Sys.time(), 'reading', opt$jx_file))
jx_info <- read.table(opt$jx_file, sep = '\t', header = FALSE,
    stringsAsFactors = FALSE, check.names = FALSE)
colnames(jx_info) <- c('chr', 'start', 'end', 'sample_ids', 'reads')

head(jx_info)


weird_df <- data.frame(chr = paste0(seqnames(weird_gr), strand(weird_gr)), start = start(weird_gr), end = end(weird_gr), stringsAsFactors = FALSE)

subset(jx_info, chr %in% weird_df$chr & start %in% weird_df$start & end %in% weird_df$end)
sapply(jx_info, class)

## Reproducibility info
proc.time()
options(width = 120)
session_info()