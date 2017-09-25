
#' linear model simulation.
#' only used on UNIX server.

linearModelSimulation<- function(n=200,
                                 p=10000,
                                 cl=20,
                                 non_zero_num_cl=5,
                                 non_zero_num_beta=10,
                                 sigma=1,
                                 overlap=1,
                                 rho=0.6,
                                 family = 'gaussian') {
  p_each_cl <- p/cl
  parameters <- list(n=n,p=p,cl=cl,p_each_cluster=p_each_cl,sigma=sigma,overlap=overlap,rho=rho,family=family)
  
  S <- .covarianceMatrix(rho,sigma,p_each_cl)
  x <- NULL
  z <- NULL
  beta.not.null <- NULL
  for(i in 1:cl){
    if(i <= non_zero_num_cl){
      #new <- sample(((i-1)*p_each_cl+1):(i*p_each_cl),10)
      #cat(new,"\n")
      new <- c(((i-1)*p_each_cl+1):((i-1)*p_each_cl+non_zero_num_beta))
      beta.not.null <- c(beta.not.null,new)
    }
    x <- cbind(x,mvrnorm(n = n,mu = rep(0,p_each_cl),Sigma = S, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
    z <- cbind(z,mvrnorm(n = n,mu = rep(0,p_each_cl),Sigma = S, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
  }
  beta_x <- rep(0, p)
  beta_z <- rep(0, p)
  
  num_non_zero <- length(beta.not.null)
  #beta_x[beta.not.null] <- 2 * ((runif(num_non_zero) > 0.5) - 0.5) * abs(runif(num_non_zero,0.5,1))
  beta_x[beta.not.null] <- abs(runif(num_non_zero,0.5,1))
  beta.not.null.z <- sample(beta.not.null,size = ceiling(overlap*num_non_zero),replace = FALSE)
  beta_z[beta.not.null.z] <- beta_x[beta.not.null.z]
  y_x <- x %*% beta_x + rnorm(n)
  y_z <- z %*% beta_z + rnorm(n)
  return(list(parameters=parameters,x=x,z=z,beta_x=beta_x,beta_z=beta_z,y_x=y_x,y_z=y_z))
}

.covarianceMatrix <- function(rho,sigma,cl){
  times <- 1:cl
  H <- abs(outer(times, times, "-"))
  R <- rho^H
  D = diag(sigma,cl);
  S = D %*% R %*% D
  return(S)
}

require(MASS)
require(foreach)
library(doParallel)
#library(doMC)
#options(cores=4)
#registerDoMC()
cl <- makeCluster(40)
registerDoParallel(cl)
sml <- 
  foreach(x=1:100) %:%
  foreach(a=seq(0,1,0.25)) %do% {
    linearModelSimulation(overlap=a)
  }
# return sml[[i]][[j]]
sml <- unlist(sml,recursive=FALSE) # flattern the list into one-level
save(sml,file="sml.rda")
stopCluster(cl)





