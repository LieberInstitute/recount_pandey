#' Custom function for creating the coverage matrix given a set of regions
#'
#' Custom function for creating the coverage matrix given a set of regions given
#' a bigwig folder and a RSE object that has the bigwig file names and the AUC 
#' information.
#'
#' @param chr A character vector with the name of the chromosome.
#' @param regions A \link[GenomicRanges]{GRanges-class} object with regions
#' for \code{chr} for which to calculate the coverage matrix.
#' @param rse
#' @param bigwig_path
#' @param chunksize A single integer vector defining the chunksize to use for
#' computing the coverage matrix. Regions will be split into different chunks
#' which can be useful when using a parallel instance as defined by 
#' \code{bpparam}.
#' @param bpparam A \link[BiocParallel]{BiocParallelParam-class} instance which
#' will be used to calculate the coverage matrix in parallel. By default, 
#' \link[BiocParallel]{SerialParam-class} will be used.
#' @param chrlen The chromosome length in base pairs. If it's \code{NULL}, the 
#' chromosome length is extracted from the Rail-RNA runs GitHub repository.
#' @param verbose If \code{TRUE} basic status updates will be printed along the 
#' way.
#' @param verboseLoad If \code{TRUE} basic status updates for loading the data
#' will be printed.

#'
#' @examples
#' ## Note that this custom function requires that derfinder 1.6.4 or newer
#' ## has been installed.
#'
#' ## Define some regions
#' library('GenomicRanges')
#' regions <- GRanges(seqnames = 'chr1', IRanges(start = c(1e4, 1e5), width = #'     1e3))
#' ## Example with Pandey's data
#' load('/dcl01/leek/data/recount_pandey/rse_gene.Rdata')
#' 
#' covMat <- custom_cov_matrix(chr = 'chr1', regions = regions, rse = rse_gene,
#' bigwig_path = '/dcl01/leek/data/sunghee_analysis/processed/coverage_bigwigs')

custom_cov_matrix <- function(chr, regions, rse, bigwig_path, chunksize = 1000,
    bpparam = NULL, chrlen = NULL, verbose = TRUE, verboseLoad = verbose) {
    ## Load required libraries
    library('SummarizedExperiment')
    library('derfinder')
    library('BiocParallel')
        
    ## Check inputs
    stopifnot(packageVersion('derfinder') >= '1.6.4')
    stopifnot(is.character(chr) & length(chr) == 1)
    stopifnot((is.numeric(chunksize) | is.integer(chunksize)) & length(chunksize) == 1)
    
    ## Windows-specific info
    if(.Platform$OS.type == 'windows') {
        warning('rtracklayer does not support importing BigWig files on Windows, so this function might not work')
    }
    
    ## Find chromosome length if absent
    if(is.null(chrlen)) {
        chrinfo <- read.table('https://raw.githubusercontent.com/nellore/runs/master/gtex/hg38.sizes',
            col.names = c('chr', 'size'), colClasses = c('character',
            'integer'))
        chrlen <- chrinfo$size[chrinfo$chr == chr]
        stopifnot(length(chrlen) == 1)
    }
    
    ## Define files to look for
    sampleFiles <- file.path(bigwig_path, colData(rse)$bigwig_file)
    names(sampleFiles) <- rownames(colData(rse))
    
    ## Define library size normalization factor
    targetSize <- 40e6 * 100
    totalMapped <- colData(rse)$auc
    mappedPerXM <- totalMapped / targetSize
    
    ## Keep only regions from the chr in question
    regions <- regions[seqnames(regions) == chr]
    
    ## Split regions into chunks
    nChunks <- length(regions) %/% chunksize
    if(length(regions) %% chunksize > 0) nChunks <- nChunks + 1
    names(regions) <- seq_len(length(regions))
    
    ## Split regions into chunks
    if(nChunks == 1) {
        regs_split <- list(regions)
        names(regs_split) <- '1'
    } else {
        regs_split <- split(regions, cut(seq_len(length(regions)),
            breaks = nChunks, labels = FALSE))
    }
    
    ## Define bpparam
    if(is.null(bpparam)) bpparam <- BiocParallel::SerialParam()
    
    ## Load coverage data
    resChunks <- lapply(regs_split, derfinder:::.railMatChrRegion,
        sampleFiles = sampleFiles, chr = chr, mappedPerXM = mappedPerXM,
        L = 1, verbose = verbose, BPPARAM.railChr = bpparam,
        verboseLoad = verboseLoad, chrlen = chrlen)
    
    ## Group results from chunks
    coverageMatrix <- do.call(rbind, lapply(resChunks, '[[', 'coverageMatrix'))
    
    ## Finish
    return(coverageMatrix)
}
