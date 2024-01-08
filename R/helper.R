#' Orthogonalized B-splines
#' @param knots Array. The knots that define the spline.
#' @param boundary_knots Array. The breakpoints that define the spline.
#' @param degree Integer. The degree of the piecewise polynomial.
#' @param predictors Array. The predictor variables with size p.
#' @param is_approx Boolean. The default is `FALSE`.
#' @return A list containing:
#'   \item{bsplines}{Matrix of orthogonalized B-splines with dimension (p, length(knots) + degree + 1)}
#'   \item{z}{Predictors used in generation}
#' @export
#' @examples
#' p <- 30
#' total_knots <- 10
#' degree <- 3
#' boundaries <- c(0, 1)
#' x <- seq(from = 0, to = 1, length.out = total_knots)
#' knots <- x[2:(total_knots - 1)]
#' predictors <- runif(p, min = 0, max = 1)
#' bsplines <- orthogonize_bspline(knots, boundaries, degree, predictors)
#' # Plot the first 5 B-splines
#' index <- order(bsplines$z)
#' original_par <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 5))
#' for (i in 1:5)
#'   plot(bsplines$z[index], bsplines$bsplines[index, i], main = i, type = "l")
#' par(original_par)
#' 
orthogonize_bspline <- function(
    knots, boundary_knots, degree, predictors = NULL, is_approx = FALSE
) {
  # Determine the number of replications
  nrep <- if (is_approx) {
    if (is.null(predictors)) 1e6 else min(1e6, length(predictors))
  } else {
    if (is.null(predictors)) 1e4 else length(predictors)
  }
  
  # Generate predictors
  x <- if (is.null(predictors)) seq(0, 1, length.out = nrep) else
    if (length(predictors) < nrep) sort(c(predictors, runif(nrep - length(predictors), min = 0, max = 1))) else
      predictors
  
  # Generate B-spline basis
  bs0 <- bs(x, knots = knots, Boundary.knots = boundary_knots, degree = degree, intercept = TRUE)
  num_basis_functions <- dim(bs0)[2]
  bsplines_matrix <- matrix(0, nrow = nrep, ncol = num_basis_functions)
  a_matrix <- diag(num_basis_functions)
  
  sqLL <- 1 / sqrt(num_basis_functions)
  bsplines_matrix[, 1] <- sqLL
  bsplines_matrix[, 2] <- sqrt(12) * sqLL * (x - 0.5)
  
  # Main loop for orthogonalization
  for (i in 3:num_basis_functions) {
    bsplines_matrix[, i] <- bs0[, i - 1]
    current_basis <- bsplines_matrix[, i]
    
    # Orthogonalization
    for (j in 1:(i - 1)) {
      a_matrix[i, j] <- -num_basis_functions * mean(current_basis * bsplines_matrix[, j])
      current_basis <- current_basis + a_matrix[i, j] * bsplines_matrix[, j]
    }
    
    # Normalize and update basis
    nee <- sqrt(mean(current_basis^2))
    bsplines_matrix[, i] <- current_basis / nee * sqLL
    a_matrix[i, ] <- a_matrix[i, ] * sqLL / nee
  }
  
  # Prepare output
  if (is.null(predictors) || !is_approx) {
    bsplines_output <- bsplines_matrix
    z_values <- x
  } else {
    idx <- sapply(predictors, function(p) which(x == p)[1])
    bsplines_output <- bsplines_matrix[idx, ]
    z_values <- predictors
  }
  
  return(list(
    bsplines = matrix(bsplines_output, ncol = length(knots) + degree + 1),
    z = z_values
  ))
}


#'
#' Internal function: Validate new locations for a qrglasso_object object
#'
#' @keywords internal
#' @param qrglasso_object An `qrglasso` class object.
#' @param top_k Integer. A matrix of the top K estimated functions.
#' @param degree Integer. Degree of the piecewise polynomial.
#' @param boundaries Array. Two boundary points.
#' @return `NULL`.
#'
check_predict_parameters <- function(qrglasso_object, top_k, degree, boundaries) {
  if (!inherits(qrglasso_object, "qrglasso")) {
    stop("Invalid object! Please enter a `qrglasso` object")
  }
  if (top_k <= 0) {
    stop("Please enter a positive top k")
  }
  if (degree <= 0) {
    stop("Please enter a positive degree")
  }
  if (length(boundaries) != 2) {
    stop("Please enter a size 2 boundaries.")
  }
  if (boundaries[1] >= boundaries[2]) {
    stop("Please input valid boundaries consisting of two elements in ascending order.")
  }
  total_knots = qrglasso_object$L - degree + 1
  if (total_knots <= 0) {
    stop("Please enter a smaller degree")
  }
}

