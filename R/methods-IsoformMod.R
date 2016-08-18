#' Main function to run isoform searching algorithm
#'
#' @param isom IsoformMod object
#' @param k Number of isoforms to start with
#' @param pilot.rep Number of random starting points to begin with
#' @importFrom matrixStats colMaxs 
#' @importFrom gtools rdirichlet 
#' @useDynLib EIA EIA_SearchMode
findIsoforms <- function(isom, k = 15, pilot.rep = 300, pilot.max.iter = 30, pilot.tol = 0.001){

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

    # pick top 5% to finish
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
