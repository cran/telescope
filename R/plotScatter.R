#' Pairwise scatter plots of the data. 
#'
#' @description Scatter plots of the observations are plotted by
#'     selecting pairs of dimensions, potentially colored by a known
#'     classification.
#' 
#' @param x A matrix or data.frame; the data consisting of metric variables.
#' @param z A vector; indicating the color to use for the observations.
#' @param label A character string; the text to include in the axes labels.
#' @param trim A scalar numeric in \eqn{[0, 0.5)}; trimming to use for quantiles to determine axes. 
#' @return `NULL`
#' @examples
#' plotScatter(iris[, 1:4], iris$Species, label = "dim")

plotScatter <- function(x, z, label = "", trim = 0) {
    stopifnot(is(x, "matrix") || is(x, "data.frame"))
    stopifnot(is(z, "vector") || is(z, "factor"))
    stopifnot(nrow(x) == length(z))
    stopifnot(is.character(label), length(label) == 1)
    stopifnot(is.numeric(trim), length(trim) == 1,
              trim >= 0, trim < 0.5)

    r <- ncol(x)
    n_plots <- r * (r-1) / 2
    if (n_plots > 25) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    cols <- pmin(5, ceiling(sqrt(n_plots)))
    rows <- pmin(5, ceiling(n_plots / cols))
    opar <- par(mfrow = c(rows, cols))
    on.exit(par(opar))

    for (j1 in 1:(r-1)) {
        for (j2 in (j1+1):(r)) {
            plot(x[, j1], x[, j2], col = z,
                 xlim = quantile(x[, j1], probs = c(trim, 1-trim)),
                 ylim = quantile(x[, j2], probs = c(trim, 1-trim)),
                 xlab = paste0(label, j1),
                 ylab = paste0(label, j2),
                 type = "p")
        }
    }
}
