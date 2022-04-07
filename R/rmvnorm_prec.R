#' Simulate from a multivariate normal 
#' 
#' @param mu a vector of means of the RV
#' @param prec the precision matrix
#' @param nsims the number of draws to take
#' @return A matrix of draws (columns) across the dimension of the RW (rows)
#' 
#' @author Aaron Osgood-Zimmerman, adapted by Taylor Okonek
#' @noRd
#' @keywords internal
rmvnorm_prec <- function(mu, prec, n.sims) {
  
  if(length(mu) == 1){
    mu <- rep(mu, nrow(prec))
  }
  
  prec <- Matrix::forceSymmetric(prec)
  
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- 1
  L <- try(Matrix::Cholesky(prec, super=TRUE), silent = TRUE)
  if (class(L) == "try-error") {
    message("Joint precision matrix non-invertible. Adding 1e-6 to diagonal")
    prec <- prec + diag(nrow(prec))*(1e-6)
    suppressWarnings(L <- Matrix::Cholesky(prec, super=FALSE))
  }
  z <- Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  x <- mu + z

  return(x)
  
}