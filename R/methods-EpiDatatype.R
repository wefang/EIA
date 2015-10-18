newEpiDatatype <- function(name, tab, has.control, grange, data.dir){
    new("EpiDatatype", name = name,
        table = tab, cells = unique(as.character(tab$cell)),
        has.control = has.control, data.dir = data.dir,
        grange = grange)
}

getSizeFactors <- function(epidt){
    if (epidt@has.control == FALSE){
        files <- epidt@table$trt.file
        total.counts <- numeric(length(files))
        for (i in 1:length(files)){
            message(paste("loading file", files[i]))
            tmp <- readBin(files[i], integer(), 1e8)
            total.counts[i] <- sum(tmp)
        }
        size.factors <- median(total.counts)/total.counts
        epidt@table$trt.sf <- size.factors
    }
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
    }
    epidt
}
