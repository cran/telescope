#' Plot multivariate categorical data. 
#' 
#' @description Plots of the level combinations of pairs of variables
#'     are created where the size of the circle indicating a level
#'     combination is proportional to the frequency of the level
#'     combination.
#' 
#' @param x A matrix or data.frame; the data consisting of categorical
#'     variables.
#' @param bg If specified, the symbols are filled with colour(s), the
#'     vector `bg` is recycled to the number of observations.  The
#'     default is to fill the symbols with grey color. 
#' @return `NULL`
#' @examples
#' with(chickwts, plotBubble(data.frame(cut(weight, 100 * 1:5), feed)))

plotBubble <- function(x, bg = "grey") {
    stopifnot(is(x, "matrix") || is(x, "data.frame"))
    r <- ncol(x)
    dim <- ceiling(sqrt(r * (r-1)/2))
    opar <- par(mfrow = c(dim, dim))
    on.exit(par(opar))
    for (j1 in 1:(r-1)) {
        for (j2 in (j1+1):r) {
            A <- as.data.frame(table(x[, j1], x[, j2]))
            A[, 1] <- as.numeric(A[, 1])
            A[, 2] <- as.numeric(A[, 2])
            A[, 3] <- A[, 3]/sum(A[, 3])       
            radius <- sqrt(A[, 3]/(4 * pi))
            symbols(A[, 1], A[, 2], circles = radius,
                    xaxt = "n", yaxt = "n", 
                    inches = FALSE, bg = bg,
                    xlab = j1, ylab = j2)
            axis(1, at = 1:max(A[, 1]), labels = 1:max(A[, 1]))
            axis(2, at = 1:max(A[, 2]), labels = 1:max(A[, 2]))
        }
    }
}

