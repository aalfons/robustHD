# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## generate random subsets
randomSubsets <- function(n, h, nsamp = 500) replicate(nsamp, sample.int(n, h))

## generate random subsets based on hyperplane
hyperplaneSubsets <- function(x, y, h, nsamp = 500) {
  d <- dim(x)
  x <- addIntercept(x)
  subsets <- randomSubsets(d[1], d[2], nsamp)
  subsets <- apply(subsets, 2, function(i, x) {
    xi <- x[i, , drop=FALSE]  # select observations
    theta <- try(hyperplane(xi), silent=TRUE)  # hyperplane coefficients or error message
    # check if point form hyperplane and add points if not
    while(inherits(theta, "try-error") && length(i) < d[1]) {
      ind <- (seq_len(d[1]))[-i]  # not yet sampled observations
      new <- if(length(ind) > 1) sample(ind, 1) else ind
      i <- c(i, new)
      xi <- rbind(xi, x[new,])  # add new point
      theta <- try(hyperplane(xi), silent=TRUE)  # hyperplane coefficients or error message
    }
    if(inherits(theta, "try-error")) {
      # singularity issues even with all data points
      sample.int(d[1], h)  # random subset
    } else {
      # compute residuals w.r.t. hyperplane
      residuals <- y - x %*% theta
      # find observations with smallest absolute residuals
      findSmallest(abs(residuals), h)
    }
  }, x)
}

## generate subsets based on lasso solutions with 3 observations
sparseSubsets <- function(x, y, lambda, h, nsamp = 500, intercept = TRUE, 
                          eps = .Machine$double.eps, use.Gram = TRUE) {
  # obtain random subsets with only 3 observations
  n <- length(y)
  subsets <- randomSubsets(n, 3, nsamp)
  # call C++ function to compute lasso fits and find observations with 
  # smallest absolute residuals
  callBackend <- getBackend()
  subsets <- callBackend("R_sparseSubsets", R_x=x, R_y=y, R_lambda=lambda, 
                         R_h=h, R_subsets=subsets, R_intercept=intercept, 
                         R_eps=eps, R_useGram=use.Gram)
}

# ## generate initial subsets
# initialSubsets <- function(x, y, h, nsamp = 500, 
#                            type = c("sparse", "hyperplane", "random"),
#                            lambda, intercept = TRUE, eps = .Machine$double.eps, 
#                            use.Gram = TRUE) {
#   ## initializations
#   d <- dim(x)
#   type <- match.arg(type)
#   if(h < d[2] && type == "hyperplane") type <- "sparse"
#   ## generate subsets
#   if(type == "sparse") {
#     # obtain random subsets with only 3 observations
#     subsets <- randomSubsets(d[1], 3, nsamp)-1
#     # call C++ function to compute lasso fits and find observations with 
#     # smallest absolute residuals
#     callBackend <- getBackend()
#     subsets <- callBackend("R_initialSubsetsSparse", R_x=x, R_y=y, 
#                            R_subsets=subsets, R_h=h, R_lambda=lambda, R_intercept=intercept, 
#                            R_eps=eps, R_useGram=use.Gram)+1
#   } else if(type == "hyperplane") {
#     ## subsets based on hyperplane
#     x <- addIntercept(x)
#     subsets <- randomSubsets(d[1], d[2], nsamp)
#     subsets <- apply(subsets, 2, function(i, x) {
#       xi <- x[i, , drop=FALSE]  # select observations
#       theta <- try(hyperplane(xi), silent=TRUE)  # hyperplane coefficients or error message
#       # check if point form hyperplane and add points if not
#       while(inherits(theta, "try-error") && length(i) < d[1]) {
#         ind <- (seq_len(d[1]))[-i]  # not yet sampled observations
#         new <- if(length(ind) > 1) sample(ind, 1) else ind
#         i <- c(i, new)
#         xi <- rbind(xi, x[new,])  # add new point
#         theta <- try(hyperplane(xi), silent=TRUE)  # hyperplane coefficients or error message
#       }
#       if(inherits(theta, "try-error")) {
#         # singularity issues even with all data points
#         sample.int(d[1], h)  # random subset
#       } else {
#         # compute residuals w.r.t. hyperplane
#         residuals <- y - x %*% theta
#         # find observations with smallest absolute residuals
#         findSmallest(abs(residuals), h)
#       }
#     }, x)
#   } else {
#     ## random subsets
#     subsets <- randomSubsets(d[1], h, nsamp)
#   }
#   ## return subsets
#   subsets
# }
