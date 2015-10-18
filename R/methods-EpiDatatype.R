newEpiDatatype <- function(name, tab, has.control,
                           input.type = c("bam", "bed", "tagAlign", "bin"),
                           bin.width, grange, data.dir){
    new("EpiDatatype", name = name,
        table = tab, cells = unique(as.character(tab$cell)),
        input.type = input.type, bin.width = bin.width,
        has.control = has.control, data.dir = data.dir,
        grange = grange)
}

convert2bin <- function(epidt, chrlen.file){
    if (epidt@has.control == FALSE){
        files <- epidt@table$trt.file
    } else {
        files  <- c(epidt@table$trt.file, epidt@table$ctrl.file)
    }
    
    # get library size and size factors when converting
    total.counts <- numeric(length(files))
    out.files <- character(length(files))
    for (i in 1:length(files)){
        message(paste("loading file", files[i]))
        out.files[i] <- file.path(epidt@data.dir, paste0(basename(files[i]), ".bin"))
        if (input.type == "tagAlign"){
            temp <- ta2bin(files[i], chrlen.file, epidt@bin.width)
        }
        if (input.type == "bam"){
            temp <- tagAlign2bin(files[i], chrlen.file, epidt@bin.width)
        }
        writeBin(as.integer(temp), out.files[i], size = 4)
        total.counts[i] <- sum(temp)
    }
    size.factors <- median(total.counts)/total.counts
    epidt@table$trt.sf <- size.factors[1:length(epidt@table$trt.file)]
    epidt@table$trt.bin <- out.files[1:length(epidt@table$trt.file)]
    if (has.control == T){
        epidt@table$ctrl.sf <- size.factors[(length(epidt@table@trt.file)+1):length(size.factors)]
        epidt@table$ctrl.bin <- out.files[(length(epidt@table@trt.file)+1):length(out.files)]
    }

    epidt
}

getSizeFactors <- function(epidt){
    if (input.type != "bin") {
        stop("only use this function when input bin counts")
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
    epidt@table$ctrl.sf <- size.factors[(length(epidt@table@trt.file)+1):length(size.factors)]

    epidt
}

generateMatrix <- function(epidt, chrlen.file, bin.width){
    datatype <- epidt@name
    tab <- epidt@table
    output.dir <- epidt@data.dir
    if (!dir.exists(output.dir)){
        dir.create(output.dir)
    }
    load.chrlen(chrlen.file, bin.width)
    message(paste("loaded chromosome lengths. \n bin width:", bin.width))

    for (chr in 1:23){
        message(paste("reading in counts for chromosome", chr))
        trt.counts.mat <- matrix(, nrow(tab), bin.counts[chr])
        for (i in 1:nrow(tab)){
            if (epidt@input.type == "bin"){
                bin.file <- tab[i, ]$trt.file
            } else {
                bin.file  <- tab[i, ]$trt.bin
            }
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
            ctrl.counts.mat <- matrix(, nrow(tab), bin.counts[chr])
            for (i in 1:nrow(tab)){
                if (epidt@input.type == "bin"){
                    bin.file <- tab[i, ]$ctrl.file
                } else {
                    bin.file  <- tab[i, ]$ctrl.bin
                }
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
