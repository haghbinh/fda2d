#' Construct a 2d-functional data object by smoothing data.
#'
#' Discrete observations on one image  are fit with a set of smooth curves, each defined by an expansion in terms of user-selected basis functions.
#' The fitting criterion is least squares.
#' @return A 2d-functional object of class 'bifd'.
#'
#' @param y Tmp1
#' @param Bx Tmp2
#' @param By Tmp3
#' @param By Tmp4
#' @importFrom fda bifd
#' @export

smooth.2db <- function(y,Bx=NULL,By=NULL,B=NULL){
  if (is.null(B)) B <- Bx%x%By
  c_vect <- solve(t(B)%*%B)%*%t(B)%*%as.vector(y)
  y_hat <- B%*%c_vect
  return(list(c_vect=c_vect, y_hat=y_hat, y_mat=matrix(y_hat,nr=nrow(y)), B=B))
}
