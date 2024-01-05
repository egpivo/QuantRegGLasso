#' Orthogonalized B-splines
#' @param knots  Array. The knots that define the spline.
#' @param boundary_knots Array. The breakpoints that define the spline.
#' @param degree Integer. The degree of the piecewise polynomial.
#' @param predictors Array. The predictor variables with size p.
#' @param is_approx Boolean. The default is `FALSE`.
#' @return his function returns a \code{list} including:
#' \itemize{
#'  \item{bsplines}{Matrix Orthogonalized B-splines with dimenstion (p, length(knots) + degree + 1}
#'  \item{z}{Predictors used in generation)}
#' }
#' 
#' @export
#' @examples
#' p <- 30
#' total_knots <- 10
#' degree <- 3
#' boundaries <- c(0, 1)
#' x <- seq(from = 0, to = 1,length.out = total_knots)
#' knots <- x[2:(total_knots - 1)]
#' predictors <- runif(p, min = 0, max = 1)
#' bsplines <- orth_bspline(knots, boundaries, degree, predictors)
#' # Plot the first 5 B-splints
#' par(mfrow=c(1, 5))
#' for(i in 1:5)
#'   plot(bsplines$z, bsplines$bsplines[i,], main=i, type="l")
#'
orth_bspline <- function(
    knots, boundary_knots, degree, predictors = NULL, is_approx = FALSE
  ) {
    if (is_approx) {
      if (is.null(predictors)) {
        nrep <- 1e6
        x <- 1:nrep
        x <- x / nrep
      } else if (length(predictors) < 1e6) {
        nrep <- 1e6
        x <-
          sort(c(predictors, runif(
            1e6 - length(predictors), min = 0, max = 1
          )))
      } else {
        nrep <- length(predictors)
        x <- predictors
      }
    } else{
      if (is.null(predictors)) {
        nrep <- 1e4
        x <- (1:nrep) / nrep
      } else{
        nrep <- length(predictors)
        x <- predictors
      }
    }
    
    bs0 <-
      bs(
        x,
        knots = knots,
        Boundary.knots = boundary_knots,
        degree = degree,
        intercept = TRUE
      )

    LL <- dim(bs0)[2]
    bs1 <- matrix(0, nrow = nrep, ncol = LL)
    a0 <- diag(LL)
    
    sqLL <- 1 / sqrt(LL)
    bs1[, 1] <- sqLL
    bs1[, 2] <- sqrt(12) * sqLL * (x - 0.5)
    
    for (ii in 3:LL)
      bs1[, ii] <- bs0[, (ii - 1)]
    
    for (jj in 3:LL) {
      ee <- as.numeric(bs1[, jj])
      
      for (ii in 1:(jj - 1))
        a0[jj, ii] <- -LL * (mean(ee * bs1[, ii]))
      
      for (ii in 1:(jj - 1))
        ee <- ee + a0[jj, ii] * as.numeric(bs1[, ii])
      
      nee <- sqrt(mean(ee ** 2))
      bs1[, jj] <- ee / nee * sqLL
      a0[jj,] <- a0[jj,] * sqLL / nee
    }
    
    if (is.null(predictors) || is_approx == FALSE) {
      bsplines <- bs1
      z <- x
    } else {
      idx <-
        sapply(1:length(predictors), function(i)
          which(x == predictors[i])[1])
      bsplines <- bs1[unlist(idx),]
      z <- predictors
    }
    
    return(
      list(
        bsplines = matrix(bsplines, ncol = length(predictors)),
        z = z
      )
    )
}
