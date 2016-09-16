## Usage:
# qrsh -l mem_free=50G,h_vmem=60G
# mkdir -p logs
# Rscript characterize_intropolis_merged.R > logs/characterize_intropolis_merged_log.txt 2>&1


## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')
library('RColorBrewer')
library('limma')

## Load data
message(paste(Sys.time(), 'loading rse_jx.Rdata'))
system.time( load('rse_jx.Rdata') )

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
n_samp_tab <- lapply(n_samples, table)
n_samp_tab
lapply(lapply(n_samples, table), function(x) { round(x / sum(x) * 100, 2) })
lapply(n_samples, summary)


## How many junctions per sample?
n_jx <- lapply(presence, colSums)
print('Number of jx per sample')
lapply(n_jx, table)
lapply(n_jx, summary)

## Annotation class for each jx
print('Exon-exon junctions by annotation class given by UCSC knownGene hg38')
lapply(rse_list, function(x) { table(rowRanges(x)$class) })
lapply(rse_list, function(x) { 
    round(table(rowRanges(x)$class) / length(rowRanges(x)) * 100, 2) })

## Maximum number coverage for the new jx
max_cov <- lapply(as.integer(names(n_samp_tab[[1]])), function(i) {
    rowMax(assays(subset(rse_new, subset = n_samples[[1]] == i), 1)$counts)
})
names(max_cov) <- as.integer(names(n_samp_tab[[1]]))
print('Maximum jx coverage for the new jx')
lapply(max_cov, summary)



## Select all jx except those from the 'annotated' class
idx <- !rowRanges(rse_present)$class %in% 'annotated'
max_cov_p <- lapply(as.integer(names(n_samp_tab[[2]])), function(i) {
    rowMax(assays(subset(rse_present, subset = n_samples[[2]] == i &idx ),
        1)$counts)
})
names(max_cov_p) <- as.integer(names(n_samp_tab[[2]]))
print('Maximum jx coverage for the jx present in Intropolis and not "annotated" by UCSC knownGene hg38')
lapply(max_cov_p, summary)

## Number of jx in UCSC knownGene hg38:
load('/dcl01/leek/data/recount-website/rse/introns_unique_with_intropolis.Rdata')
length(introns_unique)
#[1] 287062

## Are they in Pandey's data?
introns_unique$in_pandey <- countOverlaps(introns_unique, rowRanges(rse_jx), type = 'equal') > 0


venn <- vennCounts(matrix(0, ncol = 3, dimnames = list(1,
    c('Intropolis', 'knownGene', 'Pandey'))))

venn[c(3, 5, 7), ]
venn[, 4] <- c(
    0,
    length(rse_new) - sum(rowRanges(rse_new)$class == 'annotated'),
    sum(!introns_unique$in_intropolis & !introns_unique$in_pandey),
    sum(rowRanges(rse_new)$class == 'annotated'),
    81066376 - length(rse_present) - sum(introns_unique$in_intropolis & !introns_unique$in_pandey),
    length(rse_present) - sum(rowRanges(rse_present)$class == 'annotated'),
    sum(introns_unique$in_intropolis & !introns_unique$in_pandey),
    sum(rowRanges(rse_present)$class == 'annotated')
)


## Exploratory plots
my_plot <- function(info, ...) {
    plot(x = as.integer(names(info)), y = as.integer(info), ...)
}

set.seed(20160916)
pdf('exploratory_plots_merged.pdf')
## Venn diagram by Intropolis, USCS knownGene38, Pandey's data
vennDiagram(venn,
    main = "jx by Intropolis v2, UCSC knownGene hg38, Pandey's data",
    counts.col = 'blue')
## Number of samples
xx <- lapply(lapply(n_samples, table), my_plot,
    main = 'jx by presence in multiple samples', xlab = 'Number of samples',
    ylab = 'Number of jx', type = 'o', bg = 'lightblue', pch = 21)
xx <- lapply(lapply(lapply(n_samples, table), function(x) {
    round(x / sum(x) * 100, 2) }), my_plot,
    main = 'jx by presence in multiple samples', xlab = 'Number of samples',
    ylab = 'Percent of jx', ylim = c(0, 100), bg = 'lightblue', pch = 21,
    type = 'o')
## Number of jx
xx <- lapply(n_jx, function(x) {
    bar <- boxplot(x, plot = FALSE)
    boxplot(x, main = 'Number of jx by sample', col = 'lightblue',
        outline = FALSE, ylim = c(0.95, 1.05) * range(x))
    points(bar$out ~ jitter(rep(1, length(bar$out)), amount = 0.1), pch = 21)
})
## Types per sample
xx <- lapply(lapply(rse_list, function(x) { table(rowRanges(x)$class) }),
    barplot, main = 'jx by class UCSC knownGene hg38',
    col = brewer.pal(5, 'Set1'), ylab = 'Frequency')
xx <- lapply(lapply(rse_list, function(x) { 
    round(table(rowRanges(x)$class) / length(rowRanges(x)) * 100, 2) }),
    barplot, main = 'Jx by class UCSC knownGene hg38', ylim = c(0, 100),
    col = brewer.pal(5, 'Set1'), ylab = 'Percent')
dev.off()

pdf('maximum_coverage_jx_new_merged.pdf')
## Max coverage for new jx
xx <- mapply(function(x, y, z) {
    bar <- boxplot(x, plot = FALSE)
    boxplot(x, main = paste('Max cov for new jx:', y, 'samples\nn jx:', z),
        col = 'lightblue', outline = FALSE, ylim = c(0.95, 1.05) * range(x))
    points(jitter(bar$out, amount = 0.2) ~ jitter(rep(1, length(bar$out)),
        amount = 0.1), pch = 21)
    hist(x, col = 'lightblue', main = paste('New jx:', y, 'samples\nn jx:', z),
        xlab = 'Maximum coverage', breaks=seq(min(x)-0.5, max(x)+0.5, by=1))
    }, max_cov, names(max_cov), n_samp_tab[[1]])
dev.off()

pdf('maximum_coverage_jx_Intropolis_not_annotated_UCSC_merged.pdf')
## Max coverage for jx in Intropolis but not annotated by UCSC knownGene hg38
xx <- mapply(function(x, y, z) {
    bar <- boxplot(x, plot = FALSE)
    boxplot(x, main = paste(
        'jx in Intropolis but not "annotated" by knownGene:\n',
        y, 'samples; n jx:', z),
        ylab = 'Maximum coverage',
        col = 'lightblue', outline = FALSE, ylim = c(0.95, 1.05) * range(x))
    points(jitter(bar$out, amount = 0.2) ~ jitter(rep(1, length(bar$out)),
        amount = 0.1), pch = 21)
    hist(x, col = 'lightblue', main = paste(
        'jx in Intropolis but not "annotated" by knownGene:\n',
        y, 'samples; n jx:', z),
        xlab = 'Maximum coverage', breaks=seq(min(x)-0.5, max(x)+0.5, by=1))
    }, max_cov_p, names(max_cov_p), sapply(as.integer(names(n_samp_tab[[2]])),
        function(i) sum(n_samples[[2]] == i &idx) ))
dev.off()

## Reproducibility info
proc.time()
options(width = 120)
session_info()
