newEpiDatatype <- function(name, table,
                           file.type = c("bam", "bed", "tagAlign", "bin"),
                           has.control, data.dir, sig.range, chrlen.file, bin.width, chr.count){
      new("EpiDatatype", name = name,
          table = table, cells = unique(as.character(tab$cell)),
          file.type = file.type, has.control = has.control, data.dir = data.dir,
          sig.range = sig.range,
          chrlen.file = chrlen.file, bin.width = bin.width, chr.count = chr.count)
}

#' Convert alignment files to bin counts (both treatment and control), outputs are binary .bin files, currently supports tagAlign and bam files as inputs.
#' The normalizing size factors will also be calculated during the conversion.
#'
#' @export
#' @param epidt  A \code{\link{EpiDatatype}} object.
#' @importFrom tools file_path_sans_ext
#' @importFrom ghelper ta2bin bam2bin
convert2bin <- function(epidt){
      
      # process treatment and control together
      if (epidt@has.control == FALSE){
            files <- epidt@table$trt.file
      } else {
            files  <- c(epidt@table$trt.file, epidt@table$ctrl.file)
      }
      
      # get library size and size factors when converting
      total.counts <- numeric(length(files))
      out.files <- character(length(files))
      
      if (!dir.exists(epidt@data.dir)){
            dir.create(epidt@data.dir)
      }
      
      for (i in 1:length(files)){
            message(paste("converting file", files[i]))
            out.files[i] <- file.path(epidt@data.dir, paste0(sub("\\.[[:alnum:]]+$", "", basename(file_path_sans_ext(files[i]))), ".bin"))
            message(paste("output:", out.files[i]))
            if (epidt@file.type == "tagAlign"){
                  temp <- ta2bin(files[i], epidt@chrlen.file, epidt@bin.width)
            }
            if (epidt@file.type == "bam"){
                  temp <- bam2bin(files[i], epidt@chrlen.file, epidt@bin.width)
            }
            if (!epidt@file.type %in% c("bam", "tagAlign")){
                  stop(" Currently only supporting bam .bam and .tagAlign")
            }
            writeBin(as.integer(temp), out.files[i], size = 4)
            total.counts[i] <- sum(temp)
      }
      
      size.factors <- median(total.counts) / total.counts
      
      # for treatment
      # add a size factor column to the table
      epidt@table$trt.sf <- size.factors[1:length(epidt@table$trt.file)]
      # update treatment files with bin files
      epidt@table$trt.old <- epidt@table$trt.file 
      epidt@table$trt.file <- out.files[1:length(epidt@table$trt.file)]
      
      # repeat for control
      if (epidt@has.control == T){
            epidt@table$ctrl.sf <- size.factors[(length(epidt@table$trt.file)+1):length(size.factors)]
            epidt@table$ctrl.old <- epidt@table$ctrl.file
            epidt@table$ctrl.file <- out.files[(length(epidt@table$trt.file)+1):length(out.files)]
      }
      
      epidt@file.type <- "bin"
      
      epidt
}

#' Only used when file type given is binary bin counts already, otherwise, use \code{\link{convert2bin}} will calculate the size factors.
#' 
#' @param epidt  A \code{\link{EpiDatatype}} object.
getSizeFactors <- function(epidt){
      if (epidt@file.type != "bin") {
            stop("input type is not bin.")
      }
      
      if (epidt@has.control == FALSE){
            files <- epidt@table$trt.file
      } else {
            files  <- c(epidt@table$trt.file, epidt@table$ctrl.file)
      }
      
      total.counts <- numeric(length(files))
      for (i in 1:length(files)){
            message(paste("loading file", files[i]))
            tmp <- readBin(files[i], integer(), 1e8)
            total.counts[i] <- sum(tmp)
      }
      
      size.factors <- median(total.counts)/total.counts
      epidt@table$trt.sf <- size.factors[1:length(epidt@table$trt.file)]
      if (epidt@has.control == TRUE){
            epidt@table$ctrl.sf <- size.factors[(length(epidt@table$trt.file)+1):length(size.factors)]
      }
      
      epidt
}

#' Generate count matrices for each chromosome, each matrix is saved as a separate .rds files in the data directory.
#' Requires that bin counts be summarized from the alignment files by running \code{\link{convert2bin}}.
#' 
#' @export
#' @param epidt A \code{\link{EpiDatatype}} object.
#' @importFrom ghelper load.chrlen aveMatFac
#' 
generateMatrices <- function(epidt){
      
      if (epidt@file.type != "bin") {
            stop("files are not converted to .bin files yet.")
      }
      
      datatype <- epidt@name
      tab <- epidt@table
      output.dir <- epidt@data.dir
      
      if (!dir.exists(output.dir)){
            dir.create(output.dir)
      }
      
      load.chrlen(epidt@chrlen.file, epidt@bin.width)
      
      
      for (chr in 1:epidt@chr.count){
            message(paste("starting chromosome", chr, " of ", epidt$chr.count, "."))
            trt.counts.mat <- matrix(, nrow(tab), bin.counts[chr])
            for (i in 1:nrow(tab)){
                  bin.file <- tab[i, ]$trt.file
                  con <- file(bin.file, open = "rb")
                  # read only counts for current chromosome, 4 bytes integers
                  seek(con, where = 4 * bin.from[chr])
                  trt.counts.mat[i, ] <- readBin(con, integer(), bin.counts[chr])
                  close(con)
            }
            
            message("log2 transforming counts..")
            trt.conv.mat <- log2(trt.counts.mat * tab$trt.sf + 1)
            message("taking average over replicates..")
            trt.ave.mat <- aveMatFac(trt.conv.mat, tab$cell)
            message(paste("saving matrices to", output.dir))
            saveRDS(trt.counts.mat, file =
                          file.path(output.dir, paste0(datatype, "_trt_counts_chr", chr, ".rds")))
            #         saveRDS(trt.conv.mat, file =
            #                 file.path(output.dir, paste0(datatype, "_trt_conv_chr", chr, ".rds")))
            saveRDS(trt.ave.mat, file =
                          file.path(output.dir, paste0(datatype, "_trt_ave_chr", chr, ".rds")))
            
            if (epidt@has.control == TRUE){
                  message("processing control matrix...")
                  ctrl.counts.mat <- matrix(, nrow(tab), bin.counts[chr])
                  for (i in 1:nrow(tab)){
                        bin.file <- tab[i, ]$ctrl.file
                        con <- file(bin.file, open = "rb")
                        # read only counts for current chromosome, 4 bytes integers
                        seek(con, where = 4 * bin.from[chr])
                        ctrl.counts.mat[i, ] <- readBin(con, integer(), bin.counts[chr])
                        close(con)
                  }
                  message("log2 transforming counts..")
                  ctrl.conv.mat <- log2(ctrl.counts.mat * tab$ctrl.sf + 1)
                  message("taking average over replicates..")
                  ctrl.ave.mat <- aveMatFac(ctrl.conv.mat, tab$cell)
                  message(paste("saving matrices to", output.dir))
                  saveRDS(ctrl.counts.mat, file =
                                file.path(output.dir, paste0(datatype, "_ctrl_counts_chr", chr, ".rds")))
                  #             saveRDS(ctrl.conv.mat, file =
                  #                     file.path(output.dir, paste0(datatype, "_ctrl_conv_chr", chr, ".rds")))
                  saveRDS(ctrl.ave.mat, file =
                                file.path(output.dir, paste0(datatype, "_ctrl_ave_chr", chr, ".rds")))
            }
      }
      epidt
}
