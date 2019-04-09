#' Construct a 2d-functional data object by smoothing data.
#'
#' Discrete observations on one image  are fit with a set of smooth curves, each defined by an expansion in terms of user-selected basis functions.
#' The fitting criterion is least squares.
#' @return A 2d-functional object of class 'bifd'.
#'
#' @param s a set of argument values corresponding to the observations in the rows of the matrix y.
#' @param u a set of argument values corresponding to the observations in the columns of the matrix y.
#' @param y the matrix of observed values corresponding s nad u.
#' @param e a functional data basis object for the first argument s of the bivariate function.
#' @param f a functional data basis object for the second argument t of the bivariate function.
#' @examples
#' require(imager)
#' im <- as.cimg(function(x,y)
#' sin(x/16)+cos(y/16),128,128)
#' plot(im)
#' y <- im %>% as.matrix
#' m <- nrow(y)
#' n <- ncol(y)
#' p<-9
#' q<-9
#' e <- create.bspline.basis(c(0,1),p)
#' f <- create.bspline.basis(c(0,1),q)
#' u <- seq(0,1,length.out= n)
#' s <- seq(0,1,length.out= m)
#' x <- smooth.2dbasis(s,u,y,e,f)
#' y_hat <- eval.bifd(s,u,x)
#' sqrt(sum((y-y_hat)^2))
#' im_hat <- as.cimg(y_hat)
#' par(mfrow=c(1,2))
#' plot(im,main="True image")
#' plot(im_hat,main="Smoothed image")
#'
#' @importFrom fda bifd
#' @export
smooth.2dbasis <- function(s,u,y,e,f){
  p <- f$nbasis
  q <- e$nbasis
  m <- nrow(y)
  n <- ncol(y)
  E0 <- eval.basis(s,e)#m*p
  F0 <- eval.basis(u,f)#n*q
  Phi <- matrix(NA,nrow = m*n,ncol = p*q)
  for(j in 1:n) for(i in 1:m)
  {
    k <- (j-1)*m+i
    Phi[k,] <- c(E0[i,]%*%t(F0[j,]))
  }
  c_hat0 <- lsfit(Phi,c(y),intercept=FALSE)$coef
  c_mat <- matrix(c_hat0,nrow = p,ncol = q)
  # c_mat <- matrix(c(c_hat(E0,F0,c(y))),nrow = p,ncol = q)
  x <- bifd(c_mat,e,f)
  return(x)
}
