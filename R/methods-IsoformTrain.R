addDatatype <- function(isot, epidt){
    if (!is(epidt, "EpiDatatype")) stop("Please only add EpiDatatype objects")
    isot@data.types <- c(isot@data.types, epidt@name)
    isot@epidt[[epidt@name]] <- epidt
    isot
}

#' Setup the domain model
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges intersect 
#' @importFrom GenomicRanges seqnames 
#' @importFrom GenomicRanges findOverlaps 
#' @importFrom GenomicRanges subjectHits 
setupDomainModel <- function(isot, chrlen.file, bin.width, mod.dir,
                             data.types = isot@data.types){
    # domain.file is a bed file specifying the tp domain
    domain <- isot@domain
    isot@bin.width <- bin.width
    isot@mod.dir <- mod.dir
    bins.gr <- make.bins.gr(chrlen.file, bin.width)
    pb <- txtProgressBar(1, length(domain), style = 3)
    for (domain.id in 1:length(domain)){
        setTxtProgressBar(pb, domain.id)
        d <- domain[domain.id]
        chr <- as.character(seqnames(d))
        temp.mod <- new("DomainMod",
                        chr = chr,
                        start = start(d),
                        end = end(d),
                        bin.width = bin.width,
                        start.bin = bp2bin(start(d), bin.width),
                        end.bin = bp2bin(end(d), bin.width),
                        data.types = data.types)

        # find overlaped bins for each datatype
        for (j in 1:length(data.types)){
            dtype <- data.types[j]
            cells <- as.character(unique(isot@epidt[[dtype]]@table$cell))
            elements <- intersect(isot@epidt[[dtype]]@grange, d)
            # overlap of bins in the chromosome and elements
            bins <- unique(subjectHits(findOverlaps(elements,
                                             bins.gr[seqnames(bins.gr) == chr])))
            temp.mod@element.list[[dtype]] <- list(cells = cells, bins = bins)
        }
        isot@domain.list[[domain.id]] <- temp.mod
    }
    isot
}

#' generate the model files for fitting
#' @importFrom matrixStats rowSds
generateModelFile <- function(isot, id = 1:length(isot@domain)){
    if(!dir.exists(isot@mod.dir)) dir.create(isot@mod.dir)
    chr.last <- ""
    raw.mat.bydt <- list()
    if (length(id) > 1) pb <- txtProgressBar(1, length(id), style = 3)
    for (j in 1:length(id)){
        if (length(id) > 1) setTxtProgressBar(pb, j)

        # domain definition
        domain.id <- id[j]
        chr <- isot@domain.list[[domain.id]]@chr
        chr.num <- chr2num(chr)
        start.bin <- isot@domain.list[[domain.id]]@start.bin
        end.bin <- isot@domain.list[[domain.id]]@end.bin

        # for each data type make a matrix
        mat.bydt <- list()
        bg.mean.bydt <- list()
        bg.sd.bydt <- list()
        data.types <- isot@domain.list[[domain.id]]@data.types
        for (dtype in data.types) {
            bins <- isot@domain.list[[domain.id]]@element.list[[dtype]]$bins
            # prevent loading the raw matrix multiple times
            if (chr != chr.last){
                message(paste("loading data matrix for", dtype,  chr))
                raw.mat.bydt[[dtype]] <-
                    readRDS(file.path(isot@epidt[[dtype]]@data.dir,
                                      paste0(dtype, "_trt_ave_chr", chr.num, ".rds")))
            }
            if (length(bins) != 0){
                mat.bydt[[dtype]] <- raw.mat.bydt[[dtype]][, bins]
                bg.mean.bydt[[dtype]] <- matrix(rowMeans(raw.mat.bydt[[dtype]][, start.bin:end.bin]),
                                                nrow(mat.bydt[[dtype]]), length(bins))
                bg.sd.bydt[[dtype]] <- matrix(rowSds(raw.mat.bydt[[dtype]][, start.bin:end.bin]),
                                              nrow(mat.bydt[[dtype]]), length(bins))
                rownames(mat.bydt[[dtype]]) <- isot@epidt[[dtype]]@cells
                rownames(bg.mean.bydt[[dtype]]) <- isot@epidt[[dtype]]@cells
                rownames(bg.sd.bydt[[dtype]]) <- isot@epidt[[dtype]]@cells

                colnames(mat.bydt[[dtype]]) <- rep(dtype, ncol(mat.bydt[[dtype]]))
                colnames(bg.mean.bydt[[dtype]]) <- rep(dtype, ncol(mat.bydt[[dtype]])) 
                colnames(bg.sd.bydt[[dtype]]) <- rep(dtype, ncol(mat.bydt[[dtype]]))
            }
        }
        chr.last <- chr

        if (length(data.types) == 1){
            if (!is.null(mat.bydt[[data.types]])){
                mat <- mat.bydt[[data.types]]
                bg.mean <- bg.mean.bydt[[data.types]]
                bg.sd <- bg.sd.bydt[[data.types]]

                mat <- mat[complete.cases(mat), ]
                bg.mean <- bg.mean[complete.cases(bg.mean), ]
                bg.sd <- bg.sd[complete.cases(bg.sd), ]
            }
        } else {
            mat <- do.call(cbind, mat.bydt)
            bg.mean <- do.call(cbind, bg.mean.bydt)
            bg.sd <- do.call(cbind, bg.sd.bydt)
        }

        mod <- new("IsoformMod", mat = mat, bg.mean = bg.mean, bg.sd = bg.sd)
        mod.file <- file.path(isot@mod.dir, paste0("id_", domain.id, ".rda"))
        save(mod, file = mod.file)
        isot@domain.list[[domain.id]]@mod.file <- mod.file
    }
    isot
}

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

consolidateModels <- function(isot){
    pb <- txtProgressBar(1, length(isot@domain), style = 3)
    for (id in 1:length(isot@domain)){
        setTxtProgressBar(pb, id)
        mod.file <- isot@domain.list[[id]]@mod.file
        load(mod.file)
        isot@domain.list[[id]]@k <- mod@k
        isot@domain.list[[id]]@p <- mod@p
        isot@domain.list[[id]]@q <- mod@q
        isot@domain.list[[id]]@clust.like <- mod@clust.like
        isot@domain.list[[id]]@labels <- mod@labels
    }
    isot
}

fitMonocleModel <- function(isot, num.cluster, num.path, root.node = NULL, reverse = FALSE){

    clust.like.long <- do.call(cbind, lapply(1:length(isot@domain.list), function(id) {
                                                 x <- isot@domain.list[[id]]
                                                 if (dim(x@clust.like)[1] != 0){
                                                     label <- apply(x@clust.like, 1, which.max)
                                                     if (length(unique(label))!=1){
                                                         out <- x@clust.like[, sort(unique(label)), drop = F]
                                                         colnames(out) <- rep(id, ncol(out))
                                                         out
                                                     }
                                                 }
                                      }))
    #     ord2 <- hclust(dist(t(clust.like.long), method = "manhattan"))$order
    #     image.na(t(clust.like.long[, ord2]), c(0, 1))

    # temporary needs fixing
    rownames(clust.like.long) <- as.character(unique(isot@epidt[[1]]@table$cell))
    # kmeans
    kms.out <- kmeans(t(clust.like.long), num.cluster)
    centers <- t(kms.out$centers)

    d <- dist(centers)
    dist.mat <- as.matrix(d)
    ig <- graph.adjacency(dist.mat, mode = "undirected", weighted = T)
    mst <- minimum.spanning.tree(ig)

    next_node <<- 0
    res <- monocle:::pq_helper(mst, use_weights = FALSE, root_node = root.node)
    cc_ordering <- monocle:::extract_good_branched_ordering(res$subtree, res$root,
                                                                        dist.mat, num.path, reverse)
    mod <- new("MonocleMod", num.cluster = num.cluster, num.path = num.path,
                            prob.mat = clust.like.long, reduce.mat = centers, dist.mat = dist.mat,
                            mst = mst, ordering = cc_ordering)
    isot@monocle.mod <- mod
    isot
}

predictBinding <- function(isot, data.type, trt.file, peak.regions, chrlen.file, bin.width, naive = FALSE){
    load.chrlen(chrlen.file, bin.width)
    bins.gr <- make.bins.gr(chrlen.file, bin.width)
    bins.peak <- queryHits(findOverlaps(bins.gr, peak.regions))

    counts.mat <- matrix(, length(trt.file), bin.from[24])
    for (i in 1:length(trt.file)){
        counts.mat[i, ] <- readBin(trt.file[i], integer(), 3e8)
    }
    # convert counts
    total.counts <- rowSums(counts.mat)
    size.factor <- median(total.counts) / total.counts
    conv.mat <- log2(counts.mat * size.factor + 1)

    val <- colMeans(conv.mat)
    global.theta1 <- mean(val[bins.peak])
    global.sigma1 <- sd(val[bins.peak])

    output <- list()
    pb <- txtProgressBar(1, length(isot@domain), style = 3)
    for (domain.id in 1:length(isot@domain)){
        setTxtProgressBar(pb, domain.id)
        dmod <- isot@domain.list[[domain.id]]

        if (is.null(dmod@q)){
            output[[domain.id]] <- NULL
            next
        }

        chr.num <- chr2num(dmod@chr)
        bins.bydt <- lapply(dmod@element.list, function(x) x$bins)
        bins.counts.bydt <- unlist(lapply(bins.bydt, length))
        bins.from.bydt <- c(0, cumsum(bins.counts.bydt))
        dt.id <- which(names(bins.counts.bydt) == data.type)

        bins <- bins.bydt[[data.type]]
        q.sel <- sort(unique(apply(dmod@clust.like, 1, which.max)))
        tryCatch({
            q <- dmod@q[q.sel, (bins.from.bydt[[dt.id]] + 1):bins.from.bydt[[dt.id+1]]]
        }, error = function(e) {
            print(q.sel)
            print(domain.id)
            print(dim(clust.like))
        })
        if (is.null(dim(q))){
            print(domain.id)
            next
        }
        if (naive == TRUE){
            q <- matrix(rep(0.5, length(bins)), 1)
        }
        bins.gw <- chr2gw(rep(chr.num, length(bins)), bins)
        start.bin.gw <- chr2gw(chr.num, dmod@start.bin)
        end.bin.gw <- chr2gw(chr.num, dmod@end.bin)
        domain.bins.peak <- bins.peak[bins.peak > start.bin.gw & bins.peak < end.bin.gw]
        if (length(domain.bins.peak) >= 5){
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

        prob.out <- calProb(domain.val, rep(1, nrow(q))/nrow(q), q, bg.mean, bg.sd, theta1, sigma1)
        output[[domain.id]] <-
            data.frame(chr = rep(dmod@chr, length(bins)), bins = bins, prob = prob.out$a.prob)
    }
    prob.df <- do.call(rbind, output)
    gr <- with(prob.df, GRanges(seqnames = chr,
                                ranges = IRanges(start = bin2bp(bins, bin.width),
                                                 end = bin2bp(bins, bin.width) + bin.width - 1),
                                score = prob))
    return(gr)
}

#' Run differential analysis for two treatment conditions (in parallel)
#'
#' @param isot IsoformTrain object
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bplapply 
#' @importFrom matrixStats rowSds
#' @useDynLib riso riso_calLikeIso
runDifferentialAnalysis <- function(isot, data.type, bam.cond1, bam.cond2,
                                    n.perm = 2000, chrlen.file, bin.width, n.cores,
                                    id = 1:length(isot@domain)){
    multicoreParam <- MulticoreParam(workers = n.cores)

    # count matrix gw
    load.chrlen(chrlen.file, bin.width)
    n.cond1 <- length(bam.cond1)
    n.cond2 <- length(bam.cond2)

    if (n.cond1 != n.cond2) stop("don't know how to mix unequal replicates")
    counts.mat.cond1 <- matrix(, n.cond1, bin.from[24])
    counts.mat.cond2 <- matrix(, n.cond2, bin.from[24])
    for (i in 1:n.cond1){
        counts.mat.cond1[i, ] <- bam2bin(bam.cond1[i], chrlen.file, 200) 
        counts.mat.cond2[i, ] <- bam2bin(bam.cond2[i], chrlen.file, 200)
    }

    # convert matrix
    total.counts1 <- rowSums(counts.mat.cond1)
    total.counts2 <- rowSums(counts.mat.cond2)
    size.factor1 <- median(c(total.counts1, total.counts2)) / total.counts1
    size.factor2 <- median(c(total.counts1, total.counts2)) / total.counts2
    val.mat.cond1 <- log2(counts.mat.cond1 * size.factor1 + 1)
    val.mat.cond2 <- log2(counts.mat.cond2 * size.factor2 + 1)

    # average over replicates
    val.cond1 <- colMeans(val.mat.cond1)
    val.cond2 <- colMeans(val.mat.cond2)

    # estimate global theta1, sigma1
    cutoff1 <- quantile(val.cond1, 0.98)
    cutoff2 <- quantile(val.cond2, 0.98)
    theta1.global.cond1 <- mean(val.cond1[val.cond1 > cutoff1])
    theta1.global.cond2 <- mean(val.cond2[val.cond2 > cutoff2])
    sigma1.global.cond1 <- sd(val.cond1[val.cond1 > cutoff1])
    sigma1.global.cond2 <- sd(val.cond2[val.cond2 > cutoff2])

    result <-
        bplapply(1:length(id), function(i){
                     domain.id <- id[i]
                     print(domain.id)
                     dmod <- isot@domain.list[[domain.id]]
                     if (is.null(dmod@q)){
                         results[[domain.id]] <- NA
                         next
                     }

                     # subset bins and q matrix for the data type used
                     chr.num <- chr2num(dmod@chr)
                     bins.bydt <- lapply(dmod@element.list, function(x) x$bins)
                     bins.counts.bydt <- unlist(lapply(bins.bydt, length))
                     bins.from.bydt <- c(0, cumsum(bins.counts.bydt))
                     dt.id <- which(names(bins.counts.bydt) == data.type)

                     bins <- bins.bydt[[data.type]]
                     
                     # get the q matrix
                     label.tab <- table(apply(dmod@clust.like, 1, which.max))
                     q.sel <- as.numeric(names(label.tab))
                     try({
                         q <- dmod@q[q.sel, (bins.from.bydt[[dt.id]] + 1):bins.from.bydt[[dt.id+1]]]
                     })
                     if (is.null(dim(q))){
                         print(domain.id)
                         results[[domain.id]] <- NA
                         next
                     }
                     
                     tryCatch(bins.gw <- chr2gw(rep(chr.num, length(bins)), bins), warning = function(w) {print(bins.gw); print(domain.id)})
                     start.bin.gw <- chr2gw(chr.num, dmod@start.bin)
                     end.bin.gw <- chr2gw(chr.num, dmod@end.bin)

                     # setup variables to calculate likelihood for each isoform
                     domain.val.mat.cond1 <- val.mat.cond1[, bins.gw]
                     domain.val.mat.cond2 <- val.mat.cond2[, bins.gw]

                     domain.val.cond1 <- colMeans(domain.val.mat.cond1)
                     domain.val.cond2 <- colMeans(domain.val.mat.cond2)

                     theta1.est.cond1 <- theta1.global.cond1
                     theta1.est.cond2 <- theta1.global.cond2
                     sigma1.est.cond1 <- sigma1.global.cond1
                     sigma1.est.cond2 <- sigma1.global.cond2

                     #                      theta1.est.cond1 <- mean(domain.val.cond1[domain.val.cond1 > quantile(domain.val.cond1, 0.7)])
                     #                      theta1.est.cond2 <- mean(domain.val.cond2[domain.val.cond2 > quantile(domain.val.cond2, 0.7)])
                     #                      sigma1.est.cond1 <- sd(domain.val.cond1[domain.val.cond1 > quantile(domain.val.cond1, 0.7)])
                     #                      sigma1.est.cond2 <- sd(domain.val.cond2[domain.val.cond2 > quantile(domain.val.cond2, 0.7)])

                     theta1.cond1 <- rep(theta1.est.cond1, length(bins))
                     sigma1.cond1 <- rep(sigma1.est.cond1, length(bins))
                     bg.mean.cond1 <- matrix(rowMeans(val.mat.cond1[, start.bin.gw:end.bin.gw]), length(bam.cond1), length(bins))
                     bg.sd.cond1 <- matrix(rowSds(val.mat.cond1[, start.bin.gw:end.bin.gw]), length(bam.cond1), length(bins))

                     theta1.cond2 <- rep(theta1.est.cond2, length(bins))
                     sigma1.cond2 <- rep(sigma1.est.cond2, length(bins))
                     bg.mean.cond2 <- matrix(rowMeans(val.mat.cond2[, start.bin.gw:end.bin.gw]), length(bam.cond2), length(bins))
                     bg.sd.cond2 <- matrix(rowSds(val.mat.cond2[, start.bin.gw:end.bin.gw]), length(bam.cond2), length(bins))

                     # likelihood matrices
                     like.mat.cond1 <- calLikeIso(domain.val.mat.cond1, rep(1, nrow(q))/nrow(q),
                                                  q, bg.mean.cond1, bg.sd.cond1, theta1.cond1, sigma1.cond1)
                     like.mat.cond2 <- calLikeIso(domain.val.mat.cond2, rep(1, nrow(q))/nrow(q),
                                                  q, bg.mean.cond2, bg.sd.cond2, theta1.cond2, sigma1.cond2)

                     #                      max.iso <- which.max(colSums(like.mat.cond1) + colSums(like.mat.cond2))
                     #                      like.mat.cond1 <- like.mat.cond1 - like.mat.cond1[, max.iso]
                     #                      like.mat.cond2 <- like.mat.cond2 - like.mat.cond2[, max.iso]

                     like.mat.cond1 <- like.mat.cond1 - rowMeans(like.mat.cond1)
                     like.mat.cond2 <- like.mat.cond2 - rowMeans(like.mat.cond2)
                     test.stat <- colMeans(like.mat.cond1) - colMeans(like.mat.cond2)

                     # permute observations to get empirical null
                     null.stat <- matrix(, n.perm, nrow(q))
                     for(i in 1:n.perm){
                         perm.out <- permutate.mat.multi(rbind(domain.val.mat.cond1, domain.val.mat.cond2),
                                                         rbind(bg.mean.cond1, bg.mean.cond2),
                                                         rbind(bg.sd.cond1, bg.sd.cond2))
                         like.mat <- calLikeIso(perm.out[[1]], rep(1, nrow(q))/nrow(q), q, perm.out[[2]], perm.out[[3]],
                                                (theta1.cond1 + theta1.cond2)/2, (sigma1.cond1 + sigma1.cond2)/2)
                         like.mat <- like.mat - rowMeans(like.mat)
                         #                          like.mat <- like.mat - like.mat[, max.iso]
                         null.stat[i, ] <- colMeans(like.mat[1:length(bam.cond1), ]) -
                         colMeans(like.mat[(length(bam.cond1)+1):(length(bam.cond1) + length(bam.cond2)), ])
                     }

                     # calculate empirical p values
                     pvalue <- numeric(nrow(q))
                     for (k in 1:nrow(q)){
                         pvalue[k] <- empirical.pvalue(null.stat[, k], test.stat[k], "two.sided")
                     }
                     list(id = domain.id, isoform = q.sel, pvalue = pvalue, null = null.stat, test = test.stat)
}, BPPARAM = multicoreParam)

    result
}

changeModDir <- function(isot, mod.dir){
   isot@mod.dir <- mod.dir
    for (i in 1:length(isot@domain.list)){
        isot@domain.list[[i]]@mod.file <- file.path(mod.dir, basename(isot@domain.list[[i]]@mod.file))
    }
    isot
}

