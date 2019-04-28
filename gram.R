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


gram <- function(B_1 = NULL,B_2 = NULL,dx,dy, B = NULL ,method = NULL){

  if(is.null(B)==TRUE){

    B=kronecker(B_1,B_2)

  }

    if(method=="rectangular" || is.null(method)==TRUE){ # requires specification of dx and dy for integrations

      G <- dx*dy*t(B)%*%B

    } else if(method=="trapezoid"){ # future numerical integration option


    } else if(method=="simpson"){ # future numerical integration option


    }

  return(G)

}

G = gram(B = eval_b_d1, dx = delta_x, dy = 1)

# Check results

fda_pack_result <- inprod(b_d1,b_d1)

print(fda_pack_result)

print(G)

norm(G-fda_pack_result,type="F") # Frobenius norm to check accuracy

# Image of G compared to fda package results

par(mfrow=c(1,2))

image(G,main="Gram Function")

image(fda_pack_result,main="FDA package results")
