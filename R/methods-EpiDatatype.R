newEpiDatatype <- function(name, tab, has.control,
                           input.type = c("bam", "bed", "tagAlign", "bin"),
                           chrlen.file, bin.width, grange, data.dir){
    new("EpiDatatype", name = name,
        table = tab, cells = unique(as.character(tab$cell)),
        input.type = input.type, chrlen.file = chrlen.file, bin.width = bin.width,
        has.control = has.control, data.dir = data.dir,
        grange = grange)
}

convert2bin <- function(epidt){
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
        #         message(paste("loading file", files[i]))
        out.files[i] <- file.path(epidt@data.dir, paste0(sub("\\.[[:alnum:]]+$", "", basename(files[i])), ".bin"))
        #         if (epidt@input.type == "tagAlign"){
        #             temp <- ta2bin(files[i], epidt@chrlen.file, epidt@bin.width)
        #         }
        #         if (epidt@input.type == "bam"){
        #             temp <- tagAlign2bin(files[i], epidt@chrlen.file, epidt@bin.width)
        #         }
        #         if (epidt@input.type == "bin"){
        #             warning("input is already bin, conversion not needed.")
        #         }
        #         writeBin(as.integer(temp), out.files[i], size = 4)
        #         total.counts[i] <- sum(temp)
    }
    #     size.factors <- median(total.counts) / total.counts
    #     epidt@table$trt.sf <- size.factors[1:length(epidt@table$trt.file)]
    epidt@table$trt.old <- epidt@table$trt.file 
    epidt@table$trt.file <- out.files[1:length(epidt@table$trt.file)]
    if (epidt@has.control == T){
        #         epidt@table$ctrl.sf <- size.factors[(length(epidt@table$trt.file)+1):length(size.factors)]
        epidt@table$ctrl.old <- epidt@table$ctrl.file
        epidt@table$ctrl.file <- out.files[(length(epidt@table$trt.file)+1):length(out.files)]
    }

    epidt
}

getSizeFactors <- function(epidt){
    if (epidt@input.type != "bin") {
        warning("input type is not bin.")
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

generateMatrix <- function(epidt){
    datatype <- epidt@name
    tab <- epidt@table
    output.dir <- epidt@data.dir
    if (!dir.exists(output.dir)){
        dir.create(output.dir)
    }
    load.chrlen(epidt@chrlen.file, epidt@bin.width)
    message(paste("loaded chromosome lengths. \n bin width:", epidt@bin.width))

    for (chr in 1:23){
        message(paste("starting chromosome", chr))
        trt.counts.mat <- matrix(, nrow(tab), bin.counts[chr])
        for (i in 1:nrow(tab)){
            bin.file <- tab[i, ]$trt.file
            con <- file(bin.file, open = "rb")
            # read only counts for current chromosome, 4 bytes integers
            seek(con, where = 4 * bin.from[chr])
            trt.counts.mat[i, ] <- readBin(con, integer(), bin.counts[chr])
            close(con)
        }
        message("converting counts..")
        trt.conv.mat <- log2(trt.counts.mat * tab$trt.sf + 1)
        message("taking average over replicates..")
        trt.ave.mat <- aveMatFac(trt.conv.mat, tab$cell)
        message(paste("saving matrices to", output.dir))
        saveRDS(trt.counts.mat, file =
                file.path(output.dir, paste0(datatype, "_trt_counts_chr", chr, ".rds")))
        saveRDS(trt.conv.mat, file =
                file.path(output.dir, paste0(datatype, "_trt_conv_chr", chr, ".rds")))
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
            message("converting counts..")
            ctrl.conv.mat <- log2(ctrl.counts.mat * tab$ctrl.sf + 1)
            message("taking average over replicates..")
            ctrl.ave.mat <- aveMatFac(ctrl.conv.mat, tab$cell)
            message(paste("saving matrices to", output.dir))
            saveRDS(ctrl.counts.mat, file =
                    file.path(output.dir, paste0(datatype, "_ctrl_counts_chr", chr, ".rds")))
            saveRDS(ctrl.conv.mat, file =
                    file.path(output.dir, paste0(datatype, "_ctrl_conv_chr", chr, ".rds")))
            saveRDS(ctrl.ave.mat, file =
                    file.path(output.dir, paste0(datatype, "_ctrl_ave_chr", chr, ".rds")))
        }
    }
    epidt
}
