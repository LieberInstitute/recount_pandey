## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading rse_jx.Rdata'))
system.time( load('rse_jx.Rdata') )

## Get junctions in the data
message(paste(Sys.time(), 'reading jx_in_intropolis.tsv file'))
jx_intropolis_raw <- read.table('jx_in_intropolis.tsv', col.names = c('count', 'chr', 'start', 'end', 'strand'), colClasses = c('integer', 'character', 'integer', 'integer', 'character'))

## Find which are in intropolis or not
jx_intropolis <- GRanges(seqnames = jx_intropolis_raw$chr, IRanges(jx_intropolis_raw$start, jx_intropolis_raw$end), strand = jx_intropolis_raw$strand)
jx <- rowRanges(rse_jx)

jx$in_intropolis <- countOverlaps(jx, jx_intropolis, type = 'equal') > 0

print('Exon-exon junctions by presence in Intropolis')
table(jx$in_intropolis)
round(table(jx$in_intropolis) / length(jx) * 100, 2)

print('Junction class by UCSC knownGene hg38: all, new vs Intropolis, not new vs Intropolis')
table(jx$class)
table(jx$class[!jx$in_intropolis])
table(jx$class[jx$in_intropolis])

print('Junction class by UCSC knownGene hg38 (in %): all, new vs Intropolis, not new vs Intropolis')
round(table(jx$class) / length(jx$class) * 100, 2)
round(table(jx$class[!jx$in_intropolis]) / sum(!jx$in_intropolis) * 100, 2)
round(table(jx$class[jx$in_intropolis]) / sum(jx$in_intropolis) * 100, 2)

## Save results
message(paste(Sys.time(), 'saving results in rse_jx_with_intropolis.Rdata'))
rowRanges(rse_jx) <- jx
save(rse_jx, file = 'rse_jx_with_intropolis.Rdata')

## Compress original data
system('gzip jx.tsv')
system('gzip jx_in_intropolis.tsv')

## Weird ones: should be in Intropolis but are not.
weird <- jx[which(jx$class == 'annotated' & !jx$in_intropolis)]
weird

if(FALSE) {
    ## Check with Snaptron manually
    browseURL('http://stingray.cs.jhu.edu:8090/srav2/snaptron?regions=chr11:61310419-61331542&exact=1&header=1')
    browseURL('http://stingray.cs.jhu.edu:8090/srav2/snaptron?regions=chr9:70171125-70172559&exact=1&header=1')
}


## Reproducibility info
proc.time()
options(width = 120)
session_info()
