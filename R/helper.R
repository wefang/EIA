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

gatherPvalues <- function(da.out){
    do.call(rbind, lapply(da.out, function(out) {
                              if (class(out) == "try-error" | is.null(out$id)) {
                                  return(NULL)
                              } else {
                                  data.frame(id = out$id,
                                             isoform = out$isoform,
                                             pvalue = out$pvalue,
                                             statistics = out$test)
                              }
}))}
