## Usage:
# qrsh -l mem_free=50G,h_vmem=60G
# Rscript characterize_intropolis.R
# mkdir -p logs
# Rscript characterize_intropolis.R > logs/characterize_intropolis_log.txt 2>&1


## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading rse_jx_with_intropolis.Rdata'))
system.time( load('rse_jx_with_intropolis.Rdata') )

## Separate them by whether they are in Intropolis or not
rse_new <- subset(rse_jx, subset = !in_intropolis)
rse_present <- subset(rse_jx, subset = in_intropolis)
rse_list <- list('new' = rse_new, 'present' = rse_present)

## Logical matrices: jx present in sample?
presence <- lapply(rse_list, function(x) { assays(x, 1)$counts > 0 })

## Number of junctions in each case
print('Number of jx in Intropolis')
table(rowRanges(rse_jx)$in_intropolis)

print('Percent added jx vs Intropolis')
length(rowRanges(rse_new)) / 81066376 * 100

## Number of samples
n_samples <- lapply(presence, rowSums)
print('Number of samples')
lapply(n_samples, table)
lapply(lapply(n_samples, table), function(x) { round(x / sum(x) * 100, 2) })
lapply(n_samples, summary)


## How many junctions per sample?
n_jx <- lapply(presence, colSums)
print('Number of jx per sample')
lapply(n_jx, table)
lapply(n_jx, summary)

## Types per samples
print('Types per sample')
lapply(rse_list, function(x) { table(rowRanges(x)$class) })
lapply(rse_list, function(x) { round(table(rowRanges(x)$class) / length(rowRanges(x)) * 100, 2) })

## Maximum number of 
max_cov <- lapply(seq_len(10), function(i) {
    rowMax(assays(subset(rse_new, subset = n_samples[[1]] == i), 1)$counts)
})
names(max_cov) <- seq_len(10)
print('Maximum jx coverage for the new jx')
lapply(max_cov, summary)


## Exploratory plots
my_plot <- function(info, ...) {
    plot(x = as.integer(names(info)), y = as.integer(info), ...)
}

pdf('exploratory_plots.pdf')
## Number of samples
xx <- lapply(lapply(n_samples, table), my_plot, main = 'Number of samples (frequency)', xlab = 'Number of samples', ylab = 'Number of jx')
xx <- lapply(lapply(lapply(n_samples, table), function(x) { round(x / sum(x) * 100, 2) }), my_plot, main = 'Number of samples (percent)', xlab = 'Number of samples', ylab = 'Percent of jx')
## Number of jx
xx <- lapply(lapply(n_jx, table), my_plot, main = 'Number of non-zero jx (frequency)', xlab = 'Number of non-zero jx', ylab = 'Number of samples')
xx <- lapply(n_jx, boxplot, main = 'Number of non-zero jx by sample')
## Types per sample
xx <- lapply(lapply(rse_list, function(x) { table(rowRanges(x)$class) }), barplot, main = 'Jx by class (frequency)')
xx <- lapply(lapply(rse_list, function(x) { round(table(rowRanges(x)$class) / length(rowRanges(x)) * 100, 2) }), barplot, main = 'Jx by class (percent)', ylim = c(0, 100))
## Max coverage for new jx
xx <- mapply(function(x, y) boxplot(x, main = paste('Maximum cov for new jx: n-samples', y)), max_cov, names(max_cov))
dev.off()

## Reproducibility info
proc.time()
options(width = 120)
session_info()
