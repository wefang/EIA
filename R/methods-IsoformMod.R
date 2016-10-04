#' Main function to run isoform searching algorithm
#'
#' @param isom \code{\link{IsoformMod-class}} object
#' @param k Number of isoforms to start with
#' @param pilot.rep Number of random starting points to begin with
#' @importFrom matrixStats colMaxs
#' @importFrom gtools rdirichlet
#' @importFrom BiocParallel MulticoreParam bplapply
#' @useDynLib EIA EIA_SearchMode
findIsoforms <- function(isom, k, pilot.rep, pilot.max.iter, pilot.tol, max.iter, tol, n.cores = 1){
      
    multicore_param <- MulticoreParam(workers = n.cores, progressbar = T) 

    # set inital theta1, sigma1
    theta1 <- colMaxs(isom@mat) + 2 * colMaxs(isom@bg.sd)
    ind <- isom@mat > (isom@bg.mean + 2 * isom@bg.sd)
    sigma1 <- numeric(ncol(isom@mat))
    for (j in 1:ncol(isom@mat)){
        if (sum(ind[, j]) <= 2){
            sigma1[j] <- abs(theta1[j] - max(isom@mat[, j])) / 2
        } else {
            sigma1[j] <- sd(isom@mat[ind[ ,j], j])
        }
    }

    message(paste("running", k, "components..."))
    message("running pilot round...")
    param.list <- bplapply(1:pilot.rep, function(rep) {
        # random starting point for p and q
        p <- rdirichlet(1, rep(1, k))
        q <- matrix(runif(k * ncol(isom@mat)), k, ncol(isom@mat))
        temp <- SearchMode(isom@mat, isom@bg.mean, isom@bg.sd, k, pilot.max.iter, pilot.tol,
                           p, q, theta1, sigma1)
        list(loglike = temp$loglike, p = temp$p, q = temp$q, theta1 = temp$theta1, sigma1 = temp$sigma1)
    }, BPPARAM = multicore_param)

    # pick top 5% to finish
    pilot.loglike <- unlist(lapply(param.list, function(x) x$loglike))
    top <- order(pilot.loglike, decreasing = T)[1:ceiling(pilot.rep * 0.05)]
    final.param <- param.list[top]
    final.output <- bplapply(final.param, function(param) {
                               SearchMode(isom@mat, isom@bg.mean, isom@bg.sd, k, max.iter, tol,
                                          param$p, param$q, param$theta1, param$sigma1)
    }, BPPARAM = multicore_param)
    final.loglike <- unlist(lapply(final.output, function(x) x$loglike))
    out <- final.output[[which.max(final.loglike)]]

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
