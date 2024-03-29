#' @title Orthogonalized B-splines
#' @description Generate a set of orthogonalized B-splines using the Gram-Schmidt algorithm applied to the built-in function `splines::bs()`.
#' 
#' @param knots Array. The knots that define the spline.
#' @param boundary_knots Array. The breakpoints that define the spline.
#' @param degree Integer. The degree of the piecewise polynomial.
#' @param predictors Array. The predictor variables with size p.
#' @param is_approx Boolean. The default is `FALSE`.
#' @return A list containing:
#'   \item{\code{bsplines}}{Matrix of orthogonalized B-splines with dimensions \eqn{(p, \text{length}(knots) + \text{degree} + 1)}}
#'   \item{\code{z}}{Predictors used in generation}
#' @export
#' @examples
#' # Example: Generate and plot the first 5 orthogonalized B-splines
#' p <- 30
#' total_knots <- 10
#' degree <- 3
#' boundaries <- c(0, 1)
#' x <- seq(from = 0, to = 1, length.out = total_knots)
#' knots <- x[2:(total_knots - 1)]
#' predictors <- runif(p, min = 0, max = 1)
#' bsplines <- orthogonize_bspline(knots, boundaries, degree, predictors)
#' 
#' # Plot the first 5 B-splines
#' index <- order(bsplines$z)
#' original_par <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 5))
#' for (i in 1:5)
#'   plot(bsplines$z[index], bsplines$bsplines[index, i], main = i, type = "l")
#' par(original_par)
#' @export
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

#' @title Internal Function: Validate Parameters for Prediction with a `qrglasso` Object
#'
#' @description Internal function to validate parameters for predicting with a `qrglasso` class object.
#' 
#' @keywords internal
#' 
#' @param qrglasso_object A `qrglasso` class object.
#' @param metric_type Character. Metric type for gamma selection, e.g., `BIC`, `BIC-log`. Default is `BIC`.
#' @param top_k Integer. Top K estimated functions.
#' @param degree Integer. Degree of the piecewise polynomial.
#' @param boundaries Array. Two boundary points.
#' 
#' @return `NULL`.
#'
check_predict_parameters <- function(qrglasso_object, metric_type, top_k, degree, boundaries) {
  if (!inherits(qrglasso_object, "qrglasso")) {
    stop("Invalid object! Please enter a `qrglasso` object.")
  }
  if (!(metric_type %in% c("BIC", "BIC-log"))) {
    stop("Only accept types: `BIC` and `BIC-log`.")
  }
  if (top_k <= 0) {
    stop("Please enter a positive top k.")
  }
  if (degree <= 0) {
    stop("Please enter a positive degree.")
  }
  if (length(boundaries) != 2) {
    stop("Please enter a size 2 boundaries.")
  }
  if (boundaries[1] >= boundaries[2]) {
    stop("Please input valid boundaries consisting of two elements in ascending order.")
  }
  total_knots <- qrglasso_object$L - degree + 1
  if (total_knots <= 0) {
    stop("Please enter a smaller degree.")
  }
}

#' @title Internal Function: Plot Sequentially
#'
#' @description Internal function to plot ggplot2 objects sequentially.
#'
#' @keywords internal
#'
#' @param objs List. Valid ggplot2 objects to be plotted sequentially.
#'
#' @return `NULL`.
#'
plot_sequentially <- function(objs) {
  original_par <- par(no.readonly = TRUE)
  on.exit(par(original_par))
  par(ask = TRUE)
  suppressWarnings({
    for (obj in objs) {
      suppressWarnings(print(obj))
    }
  })
  par(ask = FALSE)
}

#' @title Internal Function: Plot Coefficient Function
#'
#' @description Internal function to plot coefficient functions using ggplot2.
#'
#' @keywords internal
#'
#' @param data Dataframe. A dataframe containing columns ``z``, ``coefficient``.
#' @param variate Character. A character representing the title.
#'
#' @return A ggplot object.
#' 
plot_coefficient_function <- function(data, variate) {
  default_theme <- theme_classic() +
    theme(
      text = element_text(size = 24),
      plot.title = element_text(hjust = 0.5)
    )
  result <- ggplot(data, aes(x = z, y = coefficient)) +
    geom_point(col = "#4634eb") +
    ggtitle(variate) +
    default_theme
  return(result)
}

#' @title Internal Function: Plot BIC Results w.r.t. lambda
#'
#' @description Internal function to plot BIC results with respect to lambda using ggplot2.
#'
#' @keywords internal
#'
#' @param data Dataframe. A dataframe containing columns ``lambda``, ``bic``.
#' @param variate Character. A character representing the title.
#'
#' @return A ggplot object.
#'
plot_bic_result <- function(data, variate) {
  default_theme <- theme_classic() +
    theme(
      text = element_text(size = 24),
      plot.title = element_text(hjust = 0.5)
    )
  result <- ggplot(data, aes(x = lambda, y = bic)) +
    geom_point(col = "#4634eb") +
    geom_line() +
    ggtitle(variate) +
    xlab(expression(lambda)) +
    ylab("BIC") +
    default_theme
  return(result)
}
