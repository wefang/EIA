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

addDatatype <- function(isot, epidt){
    if (!is(epidt, "EpiDatatype")) stop("please only add EpiDatatype objects")
    isot@data.types <- c(isot@data.types, epidt@name)
    isot@epidt[[epidt@name]] <- epidt
    isot
}

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
            elements <- BiocGenerics::intersect(isot@epidt[[dtype]]@grange, d)
            bins <- unique(subjectHits(findOverlaps(elements,
                                             bins.gr[seqnames(bins.gr) == chr])))
            temp.mod@element.list[[dtype]] <- list(cells = cells, bins = bins)
        }
        isot@domain.list[[domain.id]] <- temp.mod
    }
    isot
}

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
                # message(paste("loading data matrix for", dtype, chr))
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
        }

        mod <- new("IsoformMod", mat = mat, bg.mean = bg.mean, bg.sd = bg.sd)
        mod.file <- file.path(isot@mod.dir, paste0("id_", domain.id, ".rda"))
        save(mod, file = mod.file)
        isot@domain.list[[domain.id]]@mod.file <- mod.file
    }
    isot
}

findIsoforms <- function(isom, k = 15){

    mat <- isom@mat
    bg.mean <- isom@bg.mean
    bg.sd <- isom@bg.sd

    # set inital theta1, sigma1
    theta1 <- colMaxs(mat) + 2 * colMaxs(bg.sd)
    ind <- mat > (bg.mean + 2 * bg.sd)
    sigma1 <- numeric(ncol(mat))
    for (j in 1:ncol(mat)){
        if (sum(ind[, j]) <= 2){
            sigma1[j] <- abs(theta1[j] - max(mat[, j])) / 2
        } else {
            sigma1[j] <- sd(mat[ind[ ,j], j])
        }
    }

    # parameters for pilot round
    pilot.rep <- 100
    pilot.max.iter <- 30
    pilot.tol <- 0.01
    
    pilot.loglike <- numeric()
    param.list <- list()
    message(paste("running", k, "components..."))
    message("running pilot round...")
    pb <- txtProgressBar(1, pilot.rep, style = 3)
    for (rep in 1:pilot.rep){
        setTxtProgressBar(pb, rep)
        # random starting point for p and q
        p <- rdirichlet(1, rep(1, k))
        q <- matrix(runif(k * ncol(mat)), k, ncol(mat))
        temp <- SearchMode(mat, bg.mean, bg.sd, k, pilot.max.iter, pilot.tol,
                           p, q, theta1, sigma1)
        pilot.loglike[rep] <- temp$loglike
        param.list[[rep]] <- list(p = temp$p, q = temp$q, theta1 = temp$theta1, sigma1 = temp$sigma1)
    }
    # pick top 10% to finish
    top.indices <- order(pilot.loglike, decreasing = T)[1:ceiling(pilot.rep * 0.05)]
    second.param <- param.list[top.indices]
    second.loglike <- numeric()
    second.output <- list()
    for (rep in 1:length(second.param)){
        param <- second.param[[rep]]
        temp <- SearchMode(mat, bg.mean, bg.sd, k, 1000, 0.005,
                           param$p, param$q, param$theta1, param$sigma1)
        second.loglike[rep] <- temp$loglike
        second.output[[rep]] <- temp
    }
    out <- second.output[[which.max(second.loglike)]]
    
    isom@loglike <- out$loglike
    isom@labels <- apply(out$clust.like, 1, which.max)
    isom@k <- length(unique(isom@labels)) 
    isom@p <- out$p
    isom@q <- out$q
    isom@theta1 <- out$theta1
    isom@sigma1 <- out$sigma1
    isom@clust.like <- out$clust.like
    isom@cond.like <- out$cond.like
    isom
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
    }
    isot
}

fitMonocleModel <- function(isot, num.cluster, num.path, root.node = NULL, reverse = FALSE){

    clust.like.long <- do.call(cbind, lapply(isot@domain.list, function(x) {
                                                 if (dim(x@clust.like)[1] != 0){
                                                     label <- apply(x@clust.like, 1, which.max)
                                                     x@clust.like[, sort(unique(label))]
                                                 }
}))
    # temporary needs fixing
    row.names(clust.like.long) <- as.character(unique(isot@epidt[[1]]@table$cell))
    # kmeans
    kms.out <- kmeans(t(clust.like.long), num.cluster)
    centers <- t(kms.out$centers)

    d <- dist(centers)
    dist.mat <- as.matrix(d)
    ig <- graph.adjacency(dist.mat, mode = "undirected", weighted = T)
    mst <- minimum.spanning.tree(ig)

    next_node <<- 0
    res <- monocle:::pq_helper(mst, use_weights=FALSE, root_node = root.node)
    cc_ordering <- monocle:::extract_good_branched_ordering(res$subtree, res$root,
                                                            dist.mat, num.path, reverse)
    isot@monocle.mod <- new("MonocleMod", num.cluster = num.cluster, num.path = num.path,
                       prob.mat = clust.like.long, reduce.mat = centers, dist.mat = dist.mat,
                       mst = mst, ordering = cc_ordering)
    isot
}

calProb <- function(mat, p, q, bg.mean, bg.sd, theta1, sigma1){

    if(is.vector(mat)) mat <- matrix(mat, 1)
    if(is.vector(bg.mean)) bg.mean <- matrix(bg.mean, 1)
    if(is.vector(bg.sd)) bg.sd <- matrix(bg.sd, 1)

    I <- dim(mat)[1]
    J <- dim(mat)[2] 
    K <- dim(q)[1]

    like1 <- like0 <- matrix(, I, J)
    temp.like.sum <- temp.like.ratio <- array(, c(I, J, K))
    
    for (i in 1:I){
        for (j in 1:J){
            like0[i, j] <- dnorm(mat[i, j], bg.mean[i, j], bg.sd[i, j])
            like1[i, j] <- dnorm(mat[i, j], theta1[j], sigma1[j])
        }
    }
    
    for (i in 1:I){
        for (k in 1:K){
            for (j in 1:J){
                temp.like.sum[i, j, k] <- q[k, j] * like1[i, j] + (1 - q[k, j]) * like0[i, j]
                temp.like.ratio[i, j, k] <- q[k, j] * like1[i, j] /
                ( q[k, j] * like1[i, j] + (1 - q[k, j]) * like0[i, j] )
            }
        }
    }
    
    clust.like <- matrix(, I, K)
    all.like <- numeric(I)

    for (k in 1:K){
        for (i in 1:I){
            clust.like[i, k] <- log(p[k]) + sum(log(temp.like.sum[i, , k]))
        }
    }

    for (i in 1:I){
        temp <- clust.like[i, ]
        temp.max <- max(temp)
        temp <- temp - temp.max
        temp <- exp(temp)
        clust.like[i, ] <- temp / sum(temp)
        all.like[i] <- log(sum(temp)) + temp.max
    }
    loglike <- sum(all.like)

    cond.like <- array(, c(I, J, K))
    for (i in 1:I){
        for (k in 1:K){
            for (j in 1:J){
                cond.like[i, j, k] <- temp.like.ratio[i, j, k] * clust.like[i, k]
            }
        }
    }

    b.prob <- colMeans(clust.like)
    a.prob <- Reduce("+", lapply(1:K, function(k) cond.like[, , k]))

    return(list(b.prob = b.prob, a.prob = a.prob, loglike = loglike))
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

calLikeIso2 <- function(mat, p, q, bg.mean, bg.sd, theta1, sigma1){
    I <- dim(mat)[1]
    J <- dim(mat)[2] 
    K <- dim(q)[1]

    like1 <- like0 <- matrix(, I, J)
    temp.like.sum <- temp.like.ratio <- array(, c(I, J, K))

    for (i in 1:I){
        for (j in 1:J){
            like0[i, j] <- dnorm(mat[i, j], bg.mean[i, j], bg.sd[i, j])
            like1[i, j] <- dnorm(mat[i, j], theta1[j], sigma1[j])
        }
    }

    for (i in 1:I){
        for (k in 1:K){
            for (j in 1:J){
                temp.like.sum[i, j, k] <- q[k, j] * like1[i, j] + (1 - q[k, j]) * like0[i, j]
            }
        }
    }

    loglike.scale <- matrix(, I, K)
    for (k in 1:K){
        for (i in 1:I){
            loglike.scale[i, k] <- log(p[k]) + sum(log(temp.like.sum[i, , k]))
        }
    }
    for (i in 1:I){
        temp <- loglike.scale[i, ]
        loglike.scale[i, ] <- temp - min(temp)
    }
    return(loglike.scale)
}

differentialAnalysis <- function(isot, data.type, bam.cond1, bam.cond2, peaks.cond1, peaks.cond2,
                                 n.perm, mix.prop = 0, chrlen.file, bin.width, n.cores){

    multicoreParam <- MulticoreParam(workers = n.cores)

    # count matrix gw
    load.chrlen(chrlen.file, bin.width)
    n.cond1 <- length(bam.cond1)
    n.cond2 <- length(bam.cond2)
    if (mix.prop != 0){
        if (n.cond1 != n.cond2) stop("don't know how to mix unequal replicates")
        counts.mat.cond1 <- matrix(, n.cond1, bin.from[24])
        counts.mat.cond2 <- matrix(, n.cond2, bin.from[24])
        for (i in 1:n.cond1){
            align1 <- readGAlignments(bam.cond1[i])
            align2 <- readGAlignments(bam.cond2[i])
            mix.out <- mix.align(align1, align2, mix.prop, chrlen.file, bin.width)
            counts.mat.cond1[i, ] <- mix.out[[1]]
            counts.mat.cond2[i, ] <- mix.out[[2]]
        }
    } else {
        counts.mat.cond1 <- do.call(rbind, lapply(bam.cond1, function(bam) bam2bin(bam, chrlen.file, bin.width)))
        counts.mat.cond2 <- do.call(rbind, lapply(bam.cond2, function(bam) bam2bin(bam, chrlen.file, bin.width)))
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

    # find overlaping bins with peaks
    bins.gr <- make.bins.gr(chrlen.file, bin.width)
    bins.peak.cond1 <- queryHits(findOverlaps(bins.gr, peaks.cond1))
    bins.peak.cond2 <- queryHits(findOverlaps(bins.gr, peaks.cond2))

    # estimate global theta1, sigma1
    theta1.global.cond1 <- mean(val.cond1[bins.peak.cond1])
    theta1.global.cond2 <- mean(val.cond2[bins.peak.cond2])
    sigma1.global.cond1 <- sd(val.cond1[bins.peak.cond1])
    sigma1.global.cond2 <- sd(val.cond2[bins.peak.cond2])

    # estimate after mix
    theta1.global.cond12 <- mean(val.cond1[bins.peak.cond2])
    theta1.global.cond21 <- mean(val.cond2[bins.peak.cond1])
    sigma1.global.cond12 <- sd(val.cond1[bins.peak.cond2])
    sigma1.global.cond21 <- sd(val.cond2[bins.peak.cond1])

    result <- 
        bplapply(1:length(isot@domain), function(domain.id){

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
                     q.sel <- sort(unique(apply(dmod@clust.like, 1, which.max)))
                     try({
                         q <- dmod@q[q.sel, (bins.from.bydt[[dt.id]] + 1):bins.from.bydt[[dt.id+1]]]
                     })
                     if (is.null(dim(q))){
                         print(domain.id)
                         results[[domain.id]] <- NA
                         next
                     }

                     bins.gw <- chr2gw(rep(chr.num, length(bins)), bins)
                     start.bin.gw <- chr2gw(chr.num, dmod@start.bin)
                     end.bin.gw <- chr2gw(chr.num, dmod@end.bin)

                     # estimate local theta1, sigma1
                     domain.bins.peak.cond1 <- bins.peak.cond1[bins.peak.cond1 > start.bin.gw & bins.peak.cond1 < end.bin.gw]
                     domain.bins.peak.cond2 <- bins.peak.cond2[bins.peak.cond2 > start.bin.gw & bins.peak.cond2 < end.bin.gw]

                     if (length(domain.bins.peak.cond1) >= 5){
                         theta1.est.cond1 <-
                             mean(val.cond1[domain.bins.peak.cond1]) * (1 - mix.prop) + 
                             mean(val.cond1[domain.bins.peak.cond2]) * mix.prop
                         sigma1.est.cond1 <-
                             sd(val.cond1[domain.bins.peak.cond1]) * (1 - mix.prop) + 
                             sd(val.cond1[domain.bins.peak.cond2]) * mix.prop
                     } else {
                         theta1.est.cond1 <-
                             global.theta1.cond1 * (1 - mix.prop) + 
                             global.theta1.cond12 * mix.prop
                         sigma1.est.cond1 <-
                             global.sigma1.cond1 * (1 - mix.prop) +
                             global.sigma1.cond12 * mix.prop
                     }
                     # repeat for condition 2
                     if (length(domain.bins.peak.cond2) >= 5){
                         theta1.est.cond2 <-
                             mean(val.cond2[domain.bins.peak.cond2]) * (1 - mix.prop) + 
                             mean(val.cond2[domain.bins.peak.cond1]) * mix.prop
                         sigma1.est.cond2 <-
                             sd(val.cond2[domain.bins.peak.cond2]) * (1 - mix.prop) + 
                             sd(val.cond2[domain.bins.peak.cond1]) * mix.prop
                     } else {
                         theta1.est.cond2 <-
                             global.theta1.cond2 * (1 - mix.prop) + 
                             global.theta1.cond21 * mix.prop
                         sigma1.est.cond2 <-
                             global.sigma1.cond2 * (1 - mix.prop) +
                             global.sigma1.cond21 * mix.prop
                     }

                     # setup variables to calculate likelihood for each isoform
                     domain.val.mat.cond1 <- val.mat.cond1[, bins.gw]
                     domain.val.mat.cond2 <- val.mat.cond2[, bins.gw]

                     theta1.cond1 <- rep(theta1.est.cond1, length(bins))
                     sigma1.cond1 <- rep(sigma1.est.cond1, length(bins))
                     bg.mean.cond1 <- matrix(rowMeans(val.mat.cond1[, start.bin.gw:end.bin.gw]), length(bam.cond1), length(bins))
                     bg.sd.cond1 <- matrix(rowSds(val.mat.cond1[, start.bin.gw:end.bin.gw]), length(bam.cond1), length(bins))

                     theta1.cond2 <- rep(theta1.est.cond2, length(bins))
                     sigma1.cond2 <- rep(sigma1.est.cond2, length(bins))
                     bg.mean.cond2 <- matrix(rowMeans(val.mat.cond2[, start.bin.gw:end.bin.gw]), length(bam.cond2), length(bins))
                     bg.sd.cond2 <- matrix(rowSds(val.mat.cond2[, start.bin.gw:end.bin.gw]), length(bam.cond2), length(bins))

                     # relative likelihood
                     like.mat.cond1 <- calLikeIso(domain.val.mat.cond1, rep(1, nrow(q))/nrow(q),
                                                  q, bg.mean.cond1, bg.sd.cond1, theta1.cond1, sigma1.cond1)
                     like.mat.cond2 <- calLikeIso(domain.val.mat.cond2, rep(1, nrow(q))/nrow(q),
                                                  q, bg.mean.cond2, bg.sd.cond2, theta1.cond2, sigma1.cond2)
                     max.iso <- which.max(colSums(like.mat.cond1) + colSums(like.mat.cond2))
                     like.mat.cond1 <- like.mat.cond1 - like.mat.cond1[, max.iso]
                     like.mat.cond2 <- like.mat.cond2 - like.mat.cond2[, max.iso]
                     test.stat <- colMeans(like.mat.cond1) - colMeans(like.mat.cond2)

                     # permute observations to get null distribution
                     null.stat <- matrix(, n.perm, nrow(q))
                     for(i in 1:n.perm){
                         perm.out <- permutate.mat.multi(rbind(domain.val.mat.cond1, domain.val.mat.cond2),
                                                         rbind(bg.mean.cond1, bg.mean.cond2),
                                                         rbind(bg.sd.cond1, bg.sd.cond2))
                         like.mat <- calLikeIso(perm.out[[1]], rep(1, nrow(q))/nrow(q), q, perm.out[[2]], perm.out[[3]],
                                                (theta1.cond1 + theta1.cond2)/2, (sigma1.cond1 + sigma1.cond2)/2)
                         like.mat <- like.mat - like.mat[, max.iso]
                         null.stat[i, ] <- colMeans(like.mat[1:length(bam.cond1), ]) -
                         colMeans(like.mat[(length(bam.cond1)+1):(length(bam.cond1) + length(bam.cond2)), ])
                     }

                     # calculate empirical p values
                     pvalue <- numeric(nrow(q))
                     for (k in 1:nrow(q)){
                         pvalue[k] <- empirical.pvalue(null.stat[, k], test.stat[k], "two.sided")
                     }
                     list(id = domain.id, pvalue = pvalue, null = null.stat, test = test.stat)
        }, BPPARAM = multicoreParam)
    result
}

plotIsoformModel <- function(isom, zlim, main="", cell.names=""){
    # isoforms that are used
    isoform.index <- sort(unique(isom@labels))
    # subset used isoforms
    q.mat <- mod@q[isoform.index, , drop=F]
    # reorder used isoforms with hclust result
    if (nrow(q.mat) >= 2){
        q.ord <- hclust(dist(q.mat))$order  
    } else {
        q.ord <- 1 
    }
    # reorder cell labels
    label.fac <- factor(isom@labels, levels = isoform.index[q.ord])
    ord <- order(label.fac) 
    # isoform indicator	
    isoform.dum <- as.numeric(sort(label.fac))
    # find separator of isoforms
    isoform.sep <- which(diff(rev(isoform.dum)) == -1)

    # becuse color pallete only has 9 colors
    isoform.dum <- isoform.dum - floor((isoform.dum - 1) / 9) * 9

    layout(matrix(c(3, 4, 1, 2), 2), heights=c(5, 1.5), widths=c(1.3, 10))

    # original matrix
    par(mar = c(0, 1, 2, 0.8))
    image.na(t(isom@mat[ord, ]), zlim, col = brewer.pal(9, "YlGn"),
              axes=FALSE, rowsep = isoform.sep, sepcolor = "#2171b5", sepwidth = 0.05)
    title(main, line = 0.5, cex = 1)

    # isoform matrix
    par(mar = c(1, 1, 1, 0.8))
    image.na(t(q.mat[q.ord, , drop=F]), col = brewer.pal(9, "Blues"),
             zlim = c(0, 1),rowsep = 1:(length(isoform.index) - 1), axes=FALSE)
    
    # isoform indicator
    par(mar = c(0, 5, 2, 0))
    image.na(t(isoform.dum), zlim = c(1, max(isoform.dum, na.rm = T)), 
             col = brewer.pal(max(isoform.dum, na.rm = T) , "Pastel2"), rowsep = isoform.sep,
             axes = FALSE)
    mtext(text = cell.names[ord],
          side=2, line=0.3, at = seq(1, length(cell.names), by = 1), las = 1, cex = 0.5) 

    par(mar = c(1, 5, 1, 0))
    image.na(t(matrix(isom@p[isoform.index[q.ord]])), zlim = c(0, 1), col = brewer.pal(9, "PuRd"), cellnote = T, axes = F,
             rowsep = 1:(length(isoform.index) - 1))
    mtext(paste("isoform", 1:length(isoform.index)), side = 2, las = 1, at = seq(length(isoform.index), 1, by = -1), adj = 1.2, font = 2, cex = 0.8)
}

plotDomain <- function(isot, domain.id){
    d.mod <- isot@domain.list[[domain.id]]
    load(d.mod@mod.file, .GlobalEnv)
    main <- paste0("Index: ", domain.id, ", ", dim(mod@mat)[2], " bins, Region: " ,
                  d.mod@chr, ": ", d.mod@start, "-", d.mod@end)
    cells <- rownames(mod@mat)
    plotIsoformModel(mod, c(0, 8), main, cells)
    rm(mod, envir = .GlobalEnv)
}

    plotComponents <- function(isot, domain.id = 1:length(isot@domain)){
        comps <- unlist(lapply(isot@domain.list[domain.id], function(x) x@k))
        qplot(comps)
    }

    plotReducedMatrix <- function(isot){
        m <- isot@monocle.mod@reduce.mat
        ord <- order(arrange(isot@monocle.mod@ordering, sample_name)$pseudo_time)
        par(mar = c(4, 6, 4, 2))
        image.na(t(m[ord, hclust(dist(t(m)))$order]), zlim = c(min(m) - 0.1, max(m) + 0.1), axes = F)
        mtext(text = isot@monocle.mod@ordering$sample_name,
              side = 2, line = 0.3, at = seq(0, 1, length = nrow(m)), las = 1, cex = 0.6) 
    }

    plotMonocleModel <- function(isot, show.diameter = FALSE){
        mod <- isot@monocle.mod
        pca.red <- pca.reduce(mod@reduce.mat)
        pca.space.df <- data.frame(pca.red[, c(1, 2)])

        edge.list <- as.data.frame(get.edgelist(mod@mst))
        colnames(edge.list) <- c("source", "target")

        edge.df <- merge(edge.list, mod@ordering, by.x = "source", by.y = "sample_name", all = T)
        edge.df <- merge(edge.df, pca.space.df, by.x = "source", by.y = "row.names", all = T)
        edge.df <- rename(edge.df, source.PC1 = PC1, source.PC2 = PC2)
        edge.df <- merge(edge.df, pca.space.df, by.x = "target", by.y = "row.names", all = T)
        edge.df <- rename(edge.df, target.PC1 = PC1, target.PC2 = PC2)

        pca.space.with.state.df <- merge(mod@ordering, pca.space.df, by.x = "sample_name", by.y = "row.names")
        diam <- as.data.frame(as.vector(V(mod@mst)[get.diameter(mod@mst)]$name))
        colnames(diam) <- c("sample_name")
        diam <- arrange(merge(diam, pca.space.with.state.df, by = "sample_name"), pseudo_time)

        g <- ggplot(edge.df, aes(x = source.PC1, y = source.PC2))
        g <- g + geom_segment(aes(xend = target.PC1, yend= target.PC2, color = cell_state), na.rm = T)
        g <- g + geom_point(aes(color = cell_state), stat = "unique", na.rm = T)
        g <- g + geom_text(aes(label = source, color = cell_state, size = 0.4),
                           position = position_jitter(h = 0.2, w = 0.2), na.rm = T,
                           stat = "unique")
        if (show.diameter == T){
            g <- g + geom_path(aes(x = PC1, y = PC2), color = I("black"), size=0.75, data = diam) + 
            geom_point(aes(x = PC1, y = PC2, color = cell_state), size = I(1.5), data = diam) 
        }
        g
    }

