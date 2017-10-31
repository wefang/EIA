#' Constructor for IsoformTrain
#' 
#' @export
newIsoformTrain <- function(name, domain, mod.dir, chrlen.file, chr.count, bin.width){
    new("IsoformTrain", name = name, domain = domain, mod.dir = mod.dir,
        chrlen.file = chrlen.file, chr.count = chr.count, bin.width = bin.width)
}

#' Add EpiDatatype to IsoformTrain
#' 
#' @export
addEpiDatatype <- function(isot, epidt){
      if (!is(epidt, "EpiDatatype")) stop("Please only add EpiDatatype objects")
      if (epidt@bin.width != isot@bin.width | epidt@chrlen.file != isot@chrlen.file | isot@chr.count != epidt@chr.count){
            stop("Settings mismatch, make sure chromosome length files, bin widths, and chromosome numbers are the same.")
      }
      isot@data.types <- c(isot@data.types, epidt@name)
      isot@epidt[[epidt@name]] <- epidt
      isot
}

#' Setup the domain models
#' 
#' @export
#' @importFrom GenomicRanges start end intersect seqnames findOverlaps
#' @importFrom ghelper make.bins.gr bp2bin
setupDomainModel <- function(isot, data.types = isot@data.types){
      
      bins.gr <- make.bins.gr(isot@chrlen.file, isot@bin.width, isot@chr.count)
      for (domain.id in 1:length(domain)){
            message(paste("processing", domain.id, "of", length(domain)))
            d <- isot@domain[domain.id]
            temp.mod <- new("DomainMod",
                            chr = as.character(seqnames(d)),
                            start = start(d),
                            end = end(d),
                            bin.width = isot@bin.width,
                            start.bin = bp2bin(start(d), isot@bin.width),
                            end.bin = bp2bin(end(d), isot@bin.width),
                            data.types = data.types)
            
            # find signal bins for each datatype
            for (j in 1:length(data.types)){
                  dtype <- data.types[j]
                  experiments <- as.character(unique(isot@epidt[[dtype]]@table$experiment))

                  # the peak ranges in the domain
                  elements <- intersect(isot@epidt[[dtype]]@sig.range, d)
                  # the peak bins in the domain
                  bins <- unique(subjectHits(findOverlaps(elements,
                                                          bins.gr[seqnames(bins.gr) == seqnames(d)])))
                  temp.mod@element.list[[dtype]] <- list(experiments = experiments, bins = bins)
            }
            isot@domain.list[[domain.id]] <- temp.mod
      }
      isot
}

#' Generate the model files for fitting
#' 
#' @export
#' @importFrom matrixStats rowSds
#' @importFrom ghelper seq2num aveMatFac chr2gw load.chrlen
generateModelFile <- function(isot, id = 1:length(isot@domain)){

    if(!dir.exists(isot@mod.dir)) dir.create(isot@mod.dir)
    for (j in 1:length(id)){
        message(paste("processing domain", id[j]))

        load.chrlen(isot@chrlen.file, isot@bin.width)

        # domain definition
        domain.id <- id[j]
        chr <- isot@domain.list[[domain.id]]@chr
        start.bin <- isot@domain.list[[domain.id]]@start.bin
        end.bin <- isot@domain.list[[domain.id]]@end.bin

        # for each data type make matrices
        mat.l <- list()
        bg.mean.l <- list()
        bg.sd.l <- list()

        data.types <- isot@domain.list[[domain.id]]@data.types

        for (dtype in data.types) {
            if (isot@epidt[[dtype]]@has.control) stop("handling control not implementedhandling control not implemented")

            # domain count matrix
            tab <- isot@epidt[[dtype]]@table
            dmat <- matrix(, nrow(tab), (end.bin - start.bin + 1))
            for (i in 1:nrow(tab)){
                con <- file(tab[i, "trt.bin"], "rb") 
                seek(con, where = 4 * (chr2gw(chr, start.bin) - 1))
                dmat[i, ] <- readBin(con, integer(), (end.bin - start.bin + 1), size = 4)
                close(con)
            }
            dmat <- dmat * tab$trt.sf
            dmat <- log2(dmat + 1)
            dmat <- aveMatFac(dmat, tab$experiment)

            # signal bins
            bins <- isot@domain.list[[domain.id]]@element.list[[dtype]]$bins
            if (length(bins) == 0){
                warning(paste("no signals in domain", domain.id))
                write(domain.id, "error_log", append = T)
                next
            }

            mat.l[[dtype]] <- dmat[, (bins - start.bin + 1)]
            bg.mean.l[[dtype]] <- matrix(rowMeans(dmat), nrow(dmat), length(bins))
            bg.sd.l[[dtype]] <- matrix(rowSds(dmat), nrow(dmat), length(bins))
        }

        if (length(data.types) == 1){
            if (!is.null(mat.l[[data.types]])) {
                mat <- mat.l[[data.types]]
                bg.mean <- bg.mean.l[[data.types]]
                bg.sd <- bg.sd.l[[data.types]]
            }
        } else {
            stop("multiple data types to be implemented")
        }

        mod <- new("IsoformMod", mat = mat, bg.mean = bg.mean, bg.sd = bg.sd)
        mod.file <- file.path(isot@mod.dir, paste0("model_domain_", domain.id, ".rda"))
        save(mod, file = mod.file)
        isot@domain.list[[domain.id]]@mod.file <- mod.file
    }


    isot
}

#' Run model fit for a selected set of domains
#'
#' @export
#' @slot isot An \code{\link{IsoformTrain}} object.
#' @slot domain.id A vector of the indices of domains to run model fit.
runModelFit <- function(isot, domain.id, ...){
    for (id in domain.id){
        mod.file <- isot@domain.list[[id]]@mod.file
        message(paste("loading", mod.file))
        load(mod.file)
        if (dim(mod@mat)[2] > 5){
            mod <- findIsoforms(mod, ...)
        }
        save(mod, file = mod.file)
    }
}

#' Copy results from the model file to the main object.
#'
#' @export
#' @param isot An \code{\link{IsoformTrain}} object.
consolidateModels <- function(isot, id = 1:length(isot@domain)){
    for (i in id) {
        mod.file <- isot@domain.list[[i]]@mod.file
        load(mod.file)
        isot@domain.list[[i]]@k <- mod@k
        isot@domain.list[[i]]@p <- mod@p
        isot@domain.list[[i]]@q <- mod@q
        isot@domain.list[[i]]@clust.like <- mod@clust.like
        isot@domain.list[[i]]@labels <- mod@labels
    }
    isot
}

#' Predict probability of a high signal using isoform prior.
#'
#' @export
#' @importFrom ghelper load.chrlen make.bins.gr seq2num chr2gw bin2bp
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @useDynLib EIA EIA_calPostProb
predictSignal <- function(isot, id, bin.files, peaks, data.type = isot@data.types[1], neutral.prior = FALSE) {
    load.chrlen(isot@chrlen.file, isot@bin.width)
    bins.gr <- make.bins.gr(isot@chrlen.file, isot@bin.width, isot@chr.count)
    bins.peak <- queryHits(findOverlaps(bins.gr, peaks))

    # load and convert data
    counts.mat <- matrix(, length(bin.files), bin.from[24])
    for (i in 1:length(bin.files)){
        counts.mat[i, ] <- readBin(bin.files[i], integer(), 3e8)
    }
    total.counts <- rowSums(counts.mat)
    size.factor <- median(total.counts) / total.counts
    conv.mat <- log2(counts.mat * size.factor + 1)

    val <- colMeans(conv.mat)
    global.theta1 <- mean(val[bins.peak])
    global.sigma1 <- sd(val[bins.peak])

    output <- list()
    for (domain.id in id){
        message(domain.id)
        dmod <- isot@domain.list[[domain.id]]
        if (is.null(dmod@q)){
            warning(paste("no model output for domain", domain.id))
            output[[domain.id]] <- NULL
            next
        }

        chr.num <- seq2num(dmod@chr)
        bins.d <- lapply(dmod@element.list, function(x) x$bins)
        bins.counts.d <- unlist(lapply(bins.d, length))
        bins.from.d <- c(0, cumsum(bins.counts.d))
        dt.id <- which(names(bins.counts.d) == data.type)

        bins <- bins.d[[data.type]]
        if (neutral.prior) {
            q <- matrix(rep(0.5, length(bins)), 1)
        } else {
            q.sel <- sort(dmod@labels) 
            q <- dmod@q[q.sel, (bins.from.d[[dt.id]] + 1):bins.from.d[[dt.id+1]]]
        }
        bins.gw <- chr2gw(rep(chr.num, length(bins)), bins)
        start.bin.gw <- chr2gw(chr.num, dmod@start.bin)
        end.bin.gw <- chr2gw(chr.num, dmod@end.bin)
        domain.bins.peak <- bins.peak[bins.peak > start.bin.gw & bins.peak < end.bin.gw]
        if (length(domain.bins.peak) >= 5) {
            theta1.est <- mean(val[domain.bins.peak])
            sigma1.est <- sd(val[domain.bins.peak])
        } else {
            theta1.est <- global.theta1
            sigma1.est <- global.sigma1
        }
        domain.val <- val[bins.gw]

        theta1 <- rep(theta1.est, length(bins))
        sigma1 <- rep(sigma1.est, length(bins))
        bg.mean <- rep(mean(val[start.bin.gw:end.bin.gw]), length(bins))
        bg.sd <- rep(sd(val[start.bin.gw:end.bin.gw]), length(bins))

        post.out <- calPostProb(matrix(domain.val, nrow = 1), rep(1, nrow(q))/nrow(q), q,
                                matrix(bg.mean, nrow = 1), matrix(bg.sd, nrow = 1),
                                theta1, sigma1)
        output[[as.character(domain.id)]] <-
            list(chr = dmod@chr, bins = bins, a.prob = post.out$a.prob, b.prob = post.out$b.prob)
    }
    bin.prob <- do.call(rbind, lapply(output, function(x){
                                          data.frame(chr = rep(x$chr, length(x$bins)), bin_start = bin2bp(x$bins, isot@bin.width),
                                                     probability = x$a.prob[1, ])
}))
    iso.prob <- lapply(output, function(x) x$b.prob)
    list(post_bin = bin.prob, post_iso = iso.prob)
}

#' Get the matrix of relative log likelihood for domains.
#'
#' @export
#' @useDynLib EIA EIA_calLikeIso
getLikelihoodMatrix <- function(isot, id){
    like.mat.rel.list <- list()
    for (i in id) {
        load(isot@domain.list[[i]]@mod.file)
        if (dim(mod@q)[1] == 0 | is.null(mod@q)) {
            warning(paste("no results for domain", i))
            next
        }

        label.tab <- table(mod@labels)
        q.sel <- as.numeric(names(label.tab))

        like.mat <- calLikeIso(mod@mat, mod@p, mod@q[q.sel, , drop = F], mod@bg.mean, mod@bg.sd, mod@theta1, mod@sigma1)
        max.iso <- which.max(colSums(like.mat))

        like.mat.rel  <- like.mat - like.mat[, max.iso] 
        colnames(like.mat.rel) <- paste0("domain_", i, "_isoform_", q.sel)

        like.mat.rel.list[[i]] <- like.mat.rel[, -max.iso, drop = F]
    }

    like.mat.rel <- do.call(cbind, like.mat.rel.list)
    like.mat.rel.scale <- scale(like.mat.rel, center = T)
    like.mat.rel.scale
}

#' Run differential analysis for two treatment conditions (in parallel)
#'
#' @export
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bplapply 
#' @importFrom matrixStats rowSds colSds
#' @importFrom ghelper bam2bin
#' @useDynLib EIA EIA_calLikeIso
runDifferentialAnalysis <- function(isot, id, data.type, bam.cond1, bam.cond2, n.perm = 2000,  n.cores) {
    multicoreParam <- MulticoreParam(workers = n.cores)

    signal_pct_global <- 0.98
    signal_pct_domain <- 0.85
   
    # count matrix gw
    load.chrlen(isot@chrlen.file, isot@bin.width)
    n1 <- length(bam.cond1)
    n2 <- length(bam.cond2)

    if (n1 != n2) stop("don't know how to mix unequal replicates")
    counts1 <- matrix(, n1, bin.from[24])
    counts2 <- matrix(, n2, bin.from[24])
    for (i in 1:n1){
        counts1[i, ] <- bam2bin(bam.cond1[i], isot@chrlen.file, isot@chr.count, isot@bin.width) 
        counts2[i, ] <- bam2bin(bam.cond2[i], isot@chrlen.file, isot@chr.count, isot@bin.width)
    }

    # convert matrix
    total1 <- rowSums(counts1)
    total2 <- rowSums(counts2)
    sf1 <- median(c(total1, total2)) / total1
    sf2 <- median(c(total1, total2)) / total2
    val1 <- log2(counts1 * sf1 + 1)
    val2 <- log2(counts2 * sf2 + 1)

    # average over replicates
    val1ave <- colMeans(val1)
    val2ave <- colMeans(val2)

    #     thres1 <- quantile(val1ave, signal_pct_global)
    #     thres2 <- quantile(val2ave, signal_pct_global)
    #     theta1.global1 <- mean(val1ave[val1ave > thres1])
    #     theta1.global2 <- mean(val2[val2ave > thres2])
    #     sigma1.global1 <- sd(val1ave[val1ave > thres1])
    #     sigma1.global2 <- sd(val2ave[val2ave > thres2])

    result <-
        bplapply(1:length(id), function(i){
                     domain.id <- id[i]
                     message(paste0("domain ", domain.id))

                     dmod <- isot@domain.list[[domain.id]]
                     if (is.null(dmod@q)){
                         warning("no model found!")
                         return(NA)
                     }

                     chr.num <- seq2num(dmod@chr)
                     bins.d <- lapply(dmod@element.list, function(x) x$bins)
                     bins.counts.d <- unlist(lapply(bins.d, length))
                     bins.from.d <- c(0, cumsum(bins.counts.d))
                     dt.id <- which(names(bins.counts.d) == data.type)

                     bins <- bins.d[[data.type]]

                     q.sel <- sort(unique(apply(dmod@clust.like, 1, which.max)))
                     q <- dmod@q[q.sel, (bins.from.d[[dt.id]] + 1):bins.from.d[[dt.id+1]]]

                     bins.gw <- chr2gw(rep(chr.num, length(bins)), bins)
                     start.bin.gw <- chr2gw(chr.num, dmod@start.bin)
                     end.bin.gw <- chr2gw(chr.num, dmod@end.bin)

                     # setup variables to calculate likelihood for each isoform
                     val1.domain <- val1[, bins.gw]
                     val2.domain <- val2[, bins.gw]

                     val1ave.domain <- val1ave[bins.gw]
                     val2ave.domain <- val2ave[bins.gw]

                     #                      theta1.est1 <- theta1.global1
                     #                      theta1.est2 <- theta1.global2
                     #                      sigma1.est1 <- sigma1.global1
                     #                      sigma1.est2 <- sigma1.global2

                     theta1.est1 <- mean(val1ave.domain[val1ave.domain > quantile(val1ave.domain, signal_pct_domain)])
                     theta1.est2 <- mean(val2ave.domain[val2ave.domain > quantile(val2ave.domain, signal_pct_domain)])
                     sigma1.est1 <- sd(val1ave.domain[val1ave.domain > quantile(val1ave.domain, signal_pct_domain)])
                     sigma1.est2 <- sd(val2ave.domain[val2ave.domain > quantile(val2ave.domain, signal_pct_domain)])

                     theta11 <- rep(theta1.est1, length(bins))
                     sigma11 <- rep(sigma1.est1, length(bins))
                     bg.mean1 <- matrix(rowMeans(val1[, start.bin.gw:end.bin.gw]), n1, length(bins))
                     bg.sd1 <- matrix(rowSds(val1[, start.bin.gw:end.bin.gw]), n1, length(bins))

                     theta12 <- rep(theta1.est2, length(bins))
                     sigma12 <- rep(sigma1.est2, length(bins))
                     bg.mean2 <- matrix(rowMeans(val2[, start.bin.gw:end.bin.gw]), n2, length(bins))
                     bg.sd2 <- matrix(rowSds(val2[, start.bin.gw:end.bin.gw]), n2, length(bins))

                     # likelihood matrices
                     p <- rep(1, nrow(q)) / nrow(q)
                     like1 <- calLikeIso(val1.domain, p, q, bg.mean1, bg.sd1, theta11, sigma11)
                     like2 <- calLikeIso(val2.domain, p, q, bg.mean2, bg.sd2, theta12, sigma12)

                     max.iso <- which.max(colSums(like1) + colSums(like2))
                     like1 <- like1 - like1[, max.iso]
                     like2 <- like2 - like2[, max.iso]
                     test.stat <- colMeans(like1) - colMeans(like2)

                     # permute observations to get empirical null
                     null.stat <- matrix(, n.perm, nrow(q))
                     for (i in 1:n.perm) {
                         perm.out <- perm.mat.multi(rbind(val1.domain, val2.domain),
                                                    rbind(bg.mean1, bg.mean2),
                                                    rbind(bg.sd1, bg.sd2))
                         like.mat <- calLikeIso(perm.out[[1]], p, q, perm.out[[2]], perm.out[[3]],
                                                (theta11 + theta12) / 2, (sigma11 + sigma12) / 2)
                         like.mat <- like.mat - like.mat[, max.iso]
                         null.stat[i, ] <- colMeans(like.mat[1:n1, ]) - colMeans(like.mat[(n1 + 1):(n1 + n2), ])
                     }

                     # calculate empirical p values
                     pvalue <- numeric(nrow(q))
                     for (k in 1:nrow(q)){
                         pvalue[k] <- empirical.pvalue(null.stat[, k], test.stat[k], "two.sided")
                     }

                     out <- data.frame(domain = domain.id, isoform = q.sel, test_stat = test.stat,
                                       null_mean = colMeans(null.stat), null_sd = colSds(null.stat),
                                       pvalue = pvalue)
                     out[-max.iso, ]
}, BPPARAM = multicoreParam)
    do.call(rbind, result)
}

