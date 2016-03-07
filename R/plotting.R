#' Plot the isoform fitting results, called by \code{\link{plotDomain}}
#'
#' @param isom IsoformMod object
#' @param zlim Range of values in the obeserved matrix
#' @param col.sel Indices of selected columns to plot, default to all columns.
#' @param main Title of plot
#' @param cell.names Row names
#' @importFrom RColorBrewer brewer.pal
plotIsoformModel <- function(isom, zlim, col.sel = NULL, main="", cell.names=""){
    # plot only selected columns
    if (is.null(col.sel)){
        col.sel <- 1:ncol(isom@mat)
    }

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
    image.na(t(isom@mat[ord, col.sel]), zlim, col = brewer.pal(9, "YlGn"),
              axes=FALSE, rowsep = isoform.sep, sepcolor = "#2171b5", sepwidth = 0.05)
    title(main, line = 0.5, cex = 1)

    # isoform matrix
    par(mar = c(1, 1, 1, 0.8))
    image.na(t(q.mat[q.ord, col.sel, drop=F]), col = brewer.pal(9, "Blues"),
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
    layout(1)
}

plotDomain <- function(isot, domain.id, ...){
    d.mod <- isot@domain.list[[domain.id]]
    load(d.mod@mod.file, .GlobalEnv)
    main <- paste0("Index: ", domain.id, ", ", dim(mod@mat)[2], " bins, Region: " ,
                  d.mod@chr, ": ", d.mod@start, "-", d.mod@end)
    cells <- rownames(mod@mat)
    plotIsoformModel(mod, c(0, 8), main = main, cell.names = cells, ...)
    rm(mod, envir = .GlobalEnv)
}

#' Plot the distribution of number of isoforms across domains
#'
#' @importFrom ggplot2 qplot
plotComponents <- function(isot, domain.id = 1:length(isot@domain)){
    comps <- unlist(lapply(isot@domain.list[domain.id], function(x) {
                               labels.tab <- table(x@labels)
                               sum(labels.tab >= 2)
}))
        qplot(comps)
    }

plotReducedMatrix <- function(isot){
        m <- isot@monocle.mod@reduce.mat
        ord <- order(arrange(isot@monocle.mod@ordering, sample_name)$pseudo_time)
        par(mar = c(4, 6, 4, 2))
        image.na(t(m[ord, hclust(dist(t(m)))$order]), zlim = c(min(m) - 0.1, max(m) + 0.1), axes = F)
        mtext(text = isot@monocle.mod@ordering$sample_name,
              side = 2, line = 0.3, at = seq(1, length = nrow(m), by = 1), las = 1, cex = 0.6)
    }

plotMonocleModel <- function(isot, show.diameter = FALSE){
        mod <- isot@monocle.mod
        pca.red <- pca.reduce(mod@prob.mat)
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
        g <- g + geom_point(aes(color = cell_state), position = "jitter", stat = "unique", na.rm = T)
        g <- g + geom_text(aes(label = source, color = cell_state, size = 0.4),
                           position = position_jitter(h = 2, w = 2), na.rm = T,
                           stat = "unique")
        if (show.diameter == T){
            g <- g + geom_path(aes(x = PC1, y = PC2), color = I("black"), size=0.75, data = diam) +
            geom_point(aes(x = PC1, y = PC2, color = cell_state), size = I(1.5), data = diam)
        }
        g
    }

