########## Gram matrix generation techniques

# Test data
library(fda)
library(Rfssa)
d1 <- 4
d2 <- 5
n <- 100
m <- 200
min_x <- 0
max_x <- 1
min_y <- 0
max_y <- 2
delta_x <- (max_x-min_x)/n
delta_y <- (max_y-min_y)/m
delta_t <- delta_x*delta_y
range1 <- c(min_x,max_x)
range2 <- c(min_y,max_y)
b_d1 <- create.bspline.basis(range1,d1)
b_d2 <- create.bspline.basis(range2,d2)
eval_1 <- seq(0,1,length.out=n)
eval_2 <- seq(0,2,length.out=m)

eval_b_d1 <- eval.basis(eval_1,b_d1)
eval_b_d2 <- eval.basis(eval_2,b_d2)

B=kronecker(eval_b_d1,eval_b_d2)

# Gram matrix generation function

# can specify univariate bases matrices, bivariate basis matrix, delta x and delta y for integration, 
#and a tag which will specify the integration technique to be implemented


gram <- function(B_1 = NULL,B_2 = NULL,dx,dy, B = NULL ,tag = NULL){
  
  if(is.null(B)==TRUE){
    
    B=kronecker(B_1,B_2)
    
  }
  
    if(tag=="rectangular" || is.null(tag)==TRUE){ # requires specification of dx and dy for integrations
      
      G <- delta_x*delta_y*t(B)%*%B
      
    }

  return(G)

}

G = gram(B = B, dx = delta_x, dy = delta_y)

# Image of Gram matrix shows good results

image(G)

B_d1d2 <- bifd(B)

G_fda <- inprod(B)
