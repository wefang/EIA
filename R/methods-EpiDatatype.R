newEpiDatatype <- function(name, table,
                           file.type = c("bam", "bed", "tagAlign", "bin"),
                           has.control, data.dir, sig.range, chrlen.file, bin.width, chr.count){
      new("EpiDatatype", name = name,
          table = table, experiments = unique(as.character(table$experiment)),
          file.type = file.type, has.control = has.control, data.dir = data.dir,
          sig.range = sig.range,
          chrlen.file = chrlen.file, bin.width = bin.width, chr.count = chr.count)
}

#' Convert alignment files to bin counts (both treatment and control), outputs are binary .bin files, currently supports tagAlign and bam files as inputs.
#'
#' @export
#' @param epidt  A \code{\link{EpiDatatype}} object.
#' @importFrom tools file_path_sans_ext
#' @importFrom readr write_tsv
#' @importFrom ghelper ta2bin bam2bin
convert2bin <- function(epidt, sample_id = 1:nrow(epidt@table),
                        result_file = file.path(epidt@data.dir, "conversion_results.tsv")){
      
      # process treatment and control together
      if (epidt@has.control){
            files  <- c(epidt@table$trt.file[sample_id], epidt@table$ctrl.file[sample_id])
            paired <- c(epidt@table$trt.paired[sample_id], epidt@table$ctrl.paired[sample_id])
      } else {
            files <- epidt@table$trt.file[sample_id]
            paired <- epidt@table$trt.paired[sample_id]
      }
      
      # save total counts while converting
      total.counts <- numeric(length(files))
      output.files <- character(length(files))
      
      if (!dir.exists(epidt@data.dir)){
            dir.create(epidt@data.dir)
      }
      
      for (i in 1:length(files)){

            f <- files[i]
            message(paste("processing", f, ":", i, "of", length(files)))
            output.files[i] <- file.path(epidt@data.dir,
                                         paste0(sub("\\.[[:alnum:]]+$", "",
                                                    basename(file_path_sans_ext(f))), ".bin"))
            
            # download if url is given
            if (startsWith(f, "http")){
                message("downloading file...")
                download.file(f, "temp_dl", method = "wget", quiet = T)
                f <- "temp_dl"
            }

            if (epidt@file.type == "tagAlign"){
                  temp <- ta2bin(f, chrlen.file = epidt@chrlen.file,
                                  chr.count = epidt@chr.count,
                                  bin.width = epidt@bin.width)
            }

            if (epidt@file.type == "bam"){
                if (paired[i]){
                  message("sample is paired.")
                  temp <- bam2bin_paired(f, chrlen.file = epidt@chrlen.file,
                                  chr.count = epidt@chr.count,
                                  bin.width = epidt@bin.width)
                } else {
                  temp <- bam2bin(f, chrlen.file = epidt@chrlen.file,
                                  chr.count = epidt@chr.count,
                                  bin.width = epidt@bin.width)
                }
            }

            if (!epidt@file.type %in% c("bam", "tagAlign")){
                  stop(" Currently only supporting bam .bam and .tagAlign")
            }

            if (file.exists("temp_dl")) file.remove("temp_dl")
            message(paste("writing output to:", output.files[i]))
            writeBin(as.integer(temp), output.files[i], size = 4)
            total.counts[i] <- sum(temp)
      }
    
      temp_tab <- epidt@table[sample_id, ]
      temp_tab$trt.total <- total.counts[1:length(sample_id)]
      temp_tab$trt.bin <- output.files[1:length(sample_id)]
      
      if (epidt@has.control) {
          temp_tab$ctrl.total <- total.counts[(length(sample_id)+1):length(files)]
          temp_tab$ctrl.bin <- output.files[(length(sample_id)+1):length(files)]
      }

      write_tsv(temp_tab, result_file, append = T)

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
            trt.ave.mat <- aveMatFac(trt.conv.mat, tab$experiment)
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
                  ctrl.ave.mat <- aveMatFac(ctrl.conv.mat, tab$experiment)
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
