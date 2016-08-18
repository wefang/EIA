check.epidatatype <- function(object){
    errors <- character()
    if (length(object@has.control) != 1){
        msg <- paste("has.control should only have length 1")
        errors <- c(errors, msg)
        object@has.control <- object@has.control[1]
    }

    tab <- object@table
    # check if the data.frame has the specified columns
    if (object@has.control == TRUE){
        if (!all(c("cell", "trt.file", "ctrl.file") %in%
                colnames(tab))){
            msg <- "table is not in the selected format"
            errors <- c(errors, msg)
        }
    } else {
        if (!all(c("cell", "trt.file") %in%
                colnames(tab))){
            msg <- "table is not in the selected format"
            errors <- c(errors, msg)
        }
    }
    if (length(errors) == 0) TRUE else errors
}

#' An epigenomics data type to be added to the analysis.
#'
#' @slot name The name of the epigenomic data type.
#' @slot table A data frame with each row corresponding to a sample. Requires cell and trt.file columns, and ctrl.file columns if \code{has.control} is set to \code{TURE}.
#' @slot cell.types Names of the cell types used for the analysis.
#' @slot file.type The type of the alignment files, currently supporting .bam and .tagAlign.
#' @slot has.control A logical value indicating if the data type has input samples.
#' @slot data.dir The directory of data files
#' @slot sig.range A \code{\link{GenomicRanges}} object prespecifing where to look for signals in the genome.
#' @slot bin.width The length of bins to be used for summarizing counts.
#' @slot chrlen.file A text file with the length of chromosomes of the target genome.
#' @slot chr.count The number of chromosomes there is.
#' 
setClass("EpiDatatype",
         slots = c(name = "character",
                   table = "data.frame",
                   cells = "character",
                   file.type = "character",
                   has.control = "logical",
                   data.dir = "character",
                   sig.range = "GRanges",
                   bin.width = "numeric",
                   chrlen.file = "character",
                   chr.count = "numeric"),
         validity = check.epidatatype
         )

setClass("DomainMod",
         slots = c(domain.index = "numeric",
                   chr = "character",
                   start = "numeric",
                   end = "numeric",
                   bin.width = "numeric",
                   start.bin = "numeric",
                   end.bin = "numeric",
                   data.types = "character",
                   element.list = "list",
                   mod.file = "character",
                   k = "numeric",
                   p = "numeric",
                   q = "matrix",
                   clust.like = "matrix",
                   labels = "numeric")
         )

setClass("IsoformMod",
         slots = c(mat = "matrix",
                   bg.mean = "matrix",
                   bg.sd = "matrix",
                   k = "numeric",
                   p = "numeric",
                   q = "matrix",
                   theta1 = "numeric",
                   sigma1 = "numeric",
                   clust.like = "matrix",
                   cond.like = "array",
                   labels = "numeric",
                   K = "numeric",
                   BIC = "numeric", 
                   AIC = "numeric", 
                   loglike = "numeric")
         )

#' The main class for training an isoform model
#' 
#' @slot name Name of the analysis.
#' @slot data.types Name of the epigenomic data types.
#' @slot epidt A list of \code{\link{EpiDatatype}} objects.
#' @slot domain A \code{\link{GenomicRanges}} object specifying the local domains over which to individual isoform models are to be applied.
#' @slot domain.list A list of \code{\link{DomainMod}} objects, corresponding to each domains.
#' @slot mod.dir The directory for domain model files.
#' @slot cell.types Names of the cell types used for the analysis.
#' @slot bin.width The length of bins to be used for summarizing counts.
#' @slot chrlen.file A text file with the length of chromosomes of the target genome.
#' @slot chr.count The number of chromosomes there is.
#' 
setClass("IsoformTrain",
         slots = c(name = "character",
                   data.types = "character",
                   epidt = "list",
                   domain = "GRanges",
                   domain.list = "list",
                   mod.dir = "character",
                   cell.types = "character",
                   bin.width = "numeric",
                   chrlen.file = "character",
                   chr.count = "numeric")
)

setMethod("show", "IsoformTrain", function(object){
    cat("name:\n")
    print(object@name)
    cat("\n\n")
    cat("datatypes:\n")
    print(object@data.types)
    cat("\n")
    cat("domian:\n")
    print(object@domain)
    cat("\n")
    cat(paste("domain models list has", length(object@domain.list), "elements\n"))
    cat("\n")
    cat("model directory:\n")
    print(object@mod.dir)
})

