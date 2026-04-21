library(Rcpp)

sourceCpp("src/sc2col2stage_data.cpp")
sourceCpp("src/sc2col2stage_dist.cpp")

optimal_dispersion_sc2col2stage <- function(x, K = NULL) {
  #first determine if x is a dist object or a data matrix
  is_dist <- NULL
  if(inherits(x, "dist")) {
    is_dist <- TRUE
  } else if (is.matrix(x) || is.data.frame(x)) {
    is_dist <- FALSE
  } else {
    stop("Input x must be either a dist object or a data matrix.")
  }

  target_groups <- NULL
  N <- if (is_dist) attr(x, "Size") else nrow(x)
  if (is.null(K)) {

    A <- ceiling(N / 2)
    B <- floor(N / 2)
    target_groups <- c(A, B)
  } else {
    if (length(K) != 2) {
      stop("sc2col can only handle two groups, so K must be of length 2.")
    }
    if (any(!is.finite(K)) || any(K <= 0) || any(K != as.integer(K))) {
      stop("K must be a vector of positive integers.")
    }
    if (sum(K) != N) {
      stop(sprintf("The sum of the two group sizes in K must equal N=%d, but you provided sum(K)=%d.", N, sum(K)))
    }
    target_groups <- as.integer(K)
  }
  
  target_groups <- sort(target_groups, decreasing = FALSE)
  
  if (is_dist) {
    result <- solver_sc2col2stage_dist(x, N, target_groups)
  } else {
    result <- solver_sc2col2stage_data(x, N, target_groups)
  }
}





