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
        if (!all(c("cell", "rep", "trt.file", "ctrl.file") %in%
                colnames(tab))){
            msg <- "table is not in the selected format"
            errors <- c(errors, msg)
        }
    } else {
        if (!all(c("cell", "rep", "trt.file") %in%
                colnames(tab))){
            msg <- "table is not in the selected format"
            errors <- c(errors, msg)
        }
    }
    if (length(errors) == 0) TRUE else errors
}

setClass("EpiDatatype",
         slots = c(name = "character",
                   table = "data.frame",
                   cells = "character",
                   has.control = "logical",
                   data.dir = "character",
                   grange = "GRanges"),
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
                   clust.like = "matrix")
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

setClass("MonocleMod",
         slots = c(mst = "igraph",
                   num.cluster = "numeric",
                   num.path = "numeric",
                   prob.mat = "matrix",
                   reduce.mat = "matrix",
                   dist.mat = "matrix",
                   ordering = "data.frame"
                   )
         )

setClass("IsoformTrain",
         slots = c(name = "character",
                   data.types = "character",
                   bin.width = "numeric",
                   epidt = "list",
                   domain = "GRanges",
                   domain.list = "list",
                   mod.dir = "character",
                   cells = "character",
                   monocle.mod = "MonocleMod")
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

