#' Main function to run isoform searching algorithm
#'
#' @param isom \code{\link{IsoformMod-class}} object
#' @param k Number of isoforms to start with
#' @param pilot.rep Number of random starting points to begin with
#' @importFrom dplyr progress_estimated
#' @useDynLib EIA EIA_run_em
findIsoforms <- function(isom, k, pilot.rep, pilot.iter, pilot.tol, max.iter, tol, lambda, n.cores = 1){
      
    # set inital theta1, sigma1
    theta1 <- apply(isom@mat, 2, max) + 2 * apply(isom@bg.sd, 2, max)
    ind <- isom@mat > (isom@bg.mean + 2 * isom@bg.sd)
    sigma1 <- numeric(ncol(isom@mat))
    for (j in 1:ncol(isom@mat)){
        if (sum(ind[, j]) <= 2){
            sigma1[j] <- abs(theta1[j] - max(isom@mat[, j])) / 2
        } else {
            sigma1[j] <- sd(isom@mat[ind[ ,j], j])
        }
    }

    message(paste("running", k, "components."))
    message("running pilot round.")
    p <- rep(1/k, k)
    param.list <- list()
    max_cluster <- 1

    pb <- progress_estimated(pilot.rep/100)

    for (r in 1:pilot.rep){
        if (!r %% 100){
            pb$tick()$print()
        }
        # p_gamma <- rgamma(k, 1)
        # p <- p_gamma/sum(p_gamma)

        # random starting point for q
        qmat <- matrix(runif(k * ncol(isom@mat)), k, ncol(isom@mat))
        temp <- run_em(isom@mat, isom@bg.mean, isom@bg.sd,
                       k, pilot.iter, pilot.tol, n.cores,
                       p, qmat, theta1, sigma1,
                       theta1, sigma1,
                       lambda = lambda)
        num_cluster <- sum(temp$p > 0.01)
        if (num_cluster >= max_cluster) {
            param.list[[length(param.list)+1]] <- list(loglike = temp$loglike, p = temp$p, q = temp$q, theta1 = temp$theta1, sigma1 = temp$sigma1)
            max_cluster <- num_cluster
        }
    }

    # pick the ones with the most clusters to finish 
    pilot.nc <- unlist(lapply(param.list, function(x) sum(x$p > 0.01)))
    top <- which(pilot.nc == max(pilot.nc))
    message(paste0("most clusters: ", max(pilot.nc), ", with ", length(top), " to finish."))
    if (length(top) > 100) {
        pilot.loglike <- unlist(lapply(param.list[top], function(x) x$loglike))
        top <- top[order(pilot.loglike, decreasing = T)[1:100]]
    }
    final.output <- lapply(param.list[top], function(param) {
                                 run_em(isom@mat, isom@bg.mean, isom@bg.sd,
                                        k, max.iter, tol, n.cores,
                                        param$p, param$q, param$theta1, param$sigma1,
                                        theta1, sigma1,
                                        lambda = lambda)
    })
    final.loglike <- unlist(lapply(final.output, function(x) x$loglike))
    pick <- which.max(final.loglike)
    message(paste0("picked number ", pick, " in pilot."))
    out <- final.output[[pick]]

    isom@loglike <- out$loglike
    isom@labels <- apply(out$clust.like, 1, which.max)
    isom@k <- length(unique(isom@labels))
    isom@p <- out$p
    isom@q <- out$q
    isom@theta1 <- out$theta1
    isom@sigma1 <- out$sigma1
    isom@clust.like <- out$clust.like
    isom@cond.like <- out$cond.like
    if(out$converged) {
        message("reported run has converged.")
    } else {
        warning("reported run did not converge.")
    }
    isom
}
