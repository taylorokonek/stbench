#' Simulate from a multivariate normal 
#' 
#' @param mu a vector of means of the RV
#' @param prec the precision matrix
#' @param nsims the number of draws to take
#' @param constrain.idx.list a list, with each element containing the indices of the precision 
#' matrix that need to corresponding constraint in \code{A.mat.list}
#' @param A.mat.list a list, with each element containing the constraint matrix corresponding
#' to the indices in \code{constraint.idx.list}
#' @param constrain an indicator for whether to return both the constrained and unconstrained
#' draws, or only the unconstrained draws. Default is \code{TRUE}, to return both constrained
#' and unconstrained.
#' @return If \code{constrain}, a list of length two containing the unconstrained draws \code{x}
#' and constrained draws \code{x.c}. If \code{!constrain}, \code{x}.
#' 
#' @author Aaron Osgood-Zimmerman, adapted by Taylor Okonek
#' @noRd
#' @keywords internal
multiconstr_prec <- function(mu, prec, n.sims,
                             constrain.idx.list,
                             A.mat.list,
                             constrain = TRUE) {
  
  # obtain samples from joint posterior
  x <- rmvnorm_prec(mu = mu, prec = prec, n.sims = n.sims)
  
  if (!constrain) {
    print("no constraints in joint precision")
    return(x)
  } else {
    
    # create full A matrix based on constrain.idx.list and A.mat.list
    temp <- matrix(0, nrow = nrow(A.mat.list[[1]]), ncol = length(mu))
    temp[,constrain.idx.list[[1]]] <- A.mat.list[[1]]
    A <- temp
    if (length(constrain.idx.list) > 1) {
      for (i in 2:length(constrain.idx.list)) {
        temp <- matrix(0, nrow = nrow(A.mat.list[[i]]), ncol = length(mu))
        temp[,constrain.idx.list[[i]]] <- A.mat.list[[i]]
        A <- rbind(A, temp)
      }
    }
    
    e <- matrix(rep(0, nrow(A)), ncol = 1)
    
    Qinv.A <- Matrix::solve(prec, t(A))
    
    # this is how we do it for 1 draw
    app.constr <- function(x){
      x - Qinv.A %*% (solve(A %*% Qinv.A)) %*% (A %*% x - e)
    }
    
    if(n.sims == 1){
      x.c <- app.constr(x)
    }else{
      x.c <- do.call("cbind", apply(x, 2, app.constr))
    }
    
    return(list(x = x, x.c = x.c))
    
  } 
  
}