## Usage:
# qrsh -l mem_free=50G,h_vmem=60G
# mkdir -p logs
# Rscript characterize_intropolis_merged.R > logs/characterize_intropolis_merged_log.txt 2>&1


## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading rse_jx.Rdata'))
system.time( load('rse_jx.Rdata') )

## Separate them by whether they are in Intropolis or not
rse_new <- subset(rse_jx, subset = !in_intropolis)
rse_present <- subset(rse_jx, subset = in_intropolis)
rse_list <- list('new' = rse_new, 'present' = rse_present)

## Logical matrices: jx present in cell?
presence <- lapply(rse_list, function(x) { assays(x, 1)$counts > 0 })

## Number of cells
n_cells <- lapply(presence, rowSums)
print('Number of cells')
n_cell_tab <- lapply(n_cells, table)
n_cell_tab
lapply(lapply(n_cells, table), function(x) { round(x / sum(x) * 100, 2) })
lapply(n_cells, summary)


## How many junctions per cell?
n_jx <- lapply(presence, colSums)
print('Number of jx per cell')
lapply(n_jx, table)
lapply(n_jx, summary)

## Maximum number coverage for the new jx
max_cov <- lapply(as.integer(names(n_cell_tab[[1]])), function(i) {
    rowMax(assays(subset(rse_new, subset = n_cells[[1]] == i), 1)$counts)
})
names(max_cov) <- as.integer(names(n_cell_tab[[1]]))
print('Maximum jx coverage for the new jx')
lapply(max_cov, summary)

## Select all jx except those from the 'annotated' class
idx <- !rowRanges(rse_present)$class %in% 'annotated'
max_cov_p <- lapply(as.integer(names(n_cell_tab[[2]])), function(i) {
    rowMax(assays(subset(rse_present, subset = n_cells[[2]] == i &idx ),
        1)$counts)
})
names(max_cov_p) <- as.integer(names(n_cell_tab[[2]]))
print('Maximum jx coverage for the jx present in Intropolis and not "annotated" by UCSC knownGene hg38')
lapply(max_cov_p, summary)

## Exploratory plots
my_plot <- function(info, ...) {
    plot(x = as.integer(names(info)), y = as.integer(info), ...)
}

set.seed(20160916)
pdf('exploratory_plots_merged.pdf')
## Number of cells
xx <- lapply(lapply(n_cells, table), my_plot,
    main = 'jx by presence in multiple cells', xlab = 'Number of cells',
    ylab = 'Number of jx', type = 'o', bg = 'lightblue', pch = 21)
xx <- lapply(lapply(lapply(n_cells, table), function(x) {
    round(x / sum(x) * 100, 2) }), my_plot,
    main = 'jx by presence in multiple cells', xlab = 'Number of cells',
    ylab = 'Percent of jx', ylim = c(0, 100), bg = 'lightblue', pch = 21,
    type = 'o')
dev.off()

pdf('maximum_coverage_jx_new_merged.pdf')
## Max coverage for new jx
xx <- mapply(function(x, y, z) {
    bar <- boxplot(x, plot = FALSE)
    boxplot(x, main = paste('Max cov for new jx:', y, 'cells\nn jx:', z),
        col = 'lightblue', outline = FALSE, ylim = c(0.95, 1.05) * range(x))
    points(jitter(bar$out, amount = 0.2) ~ jitter(rep(1, length(bar$out)),
        amount = 0.1), pch = 21)
    hist(x, col = 'lightblue', main = paste('New jx:', y, 'cells\nn jx:', z),
        xlab = 'Maximum coverage', breaks=seq(min(x)-0.5, max(x)+0.5, by=1))
    }, max_cov, names(max_cov), n_cell_tab[[1]])
dev.off()

pdf('maximum_coverage_jx_Intropolis_not_annotated_UCSC_merged.pdf')
## Max coverage for jx in Intropolis but not annotated by UCSC knownGene hg38
xx <- mapply(function(x, y, z) {
    bar <- boxplot(x, plot = FALSE)
    boxplot(x, main = paste(
        'jx in Intropolis but not "annotated" by knownGene:\n',
        y, 'cells; n jx:', z),
        ylab = 'Maximum coverage',
        col = 'lightblue', outline = FALSE, ylim = c(0.95, 1.05) * range(x))
    points(jitter(bar$out, amount = 0.2) ~ jitter(rep(1, length(bar$out)),
        amount = 0.1), pch = 21)
    hist(x, col = 'lightblue', main = paste(
        'jx in Intropolis but not "annotated" by knownGene:\n',
        y, 'cells; n jx:', z),
        xlab = 'Maximum coverage', breaks=seq(min(x)-0.5, max(x)+0.5, by=1))
    }, max_cov_p, names(max_cov_p), sapply(as.integer(names(n_cell_tab[[2]])),
        function(i) sum(n_cells[[2]] == i &idx) ))
dev.off()

## Reproducibility info
proc.time()
options(width = 120)
session_info()
