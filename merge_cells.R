## Usage:
# qrsh -l mem_free=50G,h_vmem=60G
# module R/3.3.x
# mkdir -p logs
# Rscript merge_cells.R > logs/merge_cells_log.txt 2>&1


## Load libraries
library('SummarizedExperiment')
library('GenomicRanges')
library('devtools')

## Load data
message(paste(Sys.time(), 'loading data'))
load('rse_gene.Rdata')
load('rse_exon.Rdata')
load('rse_jx_with_intropolis.Rdata')

## Load info about cells
message(paste(Sys.time(), 'reading cell info'))
cells <- read.csv('CellMap_codes.csv', stringsAsFactors = FALSE)

## Change "No" to match syntax
cells$No <- gsub('_', '-', cells$No)


## Function for adding cell info
add_cell <- function(rse) {
    ## Match (when possible)
    map <- match(gsub('_.*', '', rownames(colData(rse))), cells$No)

    ## Add cell "No"
    colData(rse)$cell_no <- cells$No[map]
    colData(rse)$cell_no[is.na(map)] <- gsub('_.*', '', rownames(colData(rse)))[is.na(map)]

    ## Add cell code
    colData(rse)$cell_code <- cells$Code[map]
    colData(rse)$cell_code[is.na(map)] <- cells$Code[cells$Code == 'HUVEC']


    ## Add cell type
    colData(rse)$cell_type <- cells$Type.of.cell[map]
    colData(rse)$cell_type[is.na(map)] <- cells$Type.of.cell[cells$Code == 'HUVEC']
    
    ## Done
    return(rse)
}

message(paste(Sys.time(), 'Adding cell information to the rse objects'))
rse_gene <- add_cell(rse_gene)
rse_exon <- add_cell(rse_exon)
rse_jx <- add_cell(rse_jx)

message(paste(Sys.time(), 'Saving rse objects with cell info at rse_with_cell'))
dir.create('rse_with_cell', showWarnings = FALSE)
save(rse_gene, file = 'rse_with_cell/rse_gene.Rdata')
save(rse_exon, file = 'rse_with_cell/rse_exon.Rdata')
save(rse_jx, file = 'rse_with_cell/rse_jx.Rdata')


merge_by_cell <- function(rse) {
    ## Split by cell type
    cells_no <- unique(colData(rse)$cell_no)
    rse_split <- lapply(cells_no, function(cell) {
        ## Select corresponding technical replicates
        current <- subset(rse, select = colData(rse)$cell_no == cell)
        
        ## Build new counts by adding across replicates
        counts <- matrix(rowSums(assays(current, 1)$counts), ncol = 1)
        colnames(counts) <- cell
        
        ## Build new metadata info
        meta <- colData(current)
        metadata <- meta[1, ]
        metadata$reads_downloaded <- sum(meta$reads_downloaded)
        metadata$mapped_read_count <- sum(meta$mapped_read_count)
        metadata$auc <- sum(meta$auc)
        metadata$bigwig_file <- CharacterList(meta$bigwig_file)
        rownames(metadata) <- cell
        
        ## Finish by constructing the new object
        res <- SummarizedExperiment(assays = list('counts' = counts),
            colData = metadata, rowRanges = rowRanges(current))
            
        return(res)
    })
    names(rse_split) <- cells_no
    
    ## Combine resulting objects
    result <- do.call(cbind, rse_split)
    
    ## Done
    return(result)
}

message(paste(Sys.time(), 'Creating new rse objects with replicates merged'))
rse_gene <- merge_by_cell(rse_gene)
rse_exon <- merge_by_cell(rse_exon)
rse_jx <- merge_by_cell(rse_jx)

message(paste(Sys.time(), 'Saving rse objects with at rse_merged'))
dir.create('rse_merged', showWarnings = FALSE)
save(rse_gene, file = 'rse_merged/rse_gene.Rdata')
save(rse_exon, file = 'rse_merged/rse_exon.Rdata')
save(rse_jx, file = 'rse_merged/rse_jx.Rdata')


## Reproducibility info
proc.time()
options(width = 120)
session_info()
