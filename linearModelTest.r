setwd("~/Projects/UIC/SKI/Rcode")
source("SKI.r")
require(MASS)
require(glmnet)
require(SIS)
source("screening.r")

#' A test function for linear models

.covarianceMatrix <- function(rho,sigma,cl){
  times <- 1:cl
  H <- abs(outer(times, times, "-"))
  R <- rho^H
  D = diag(sigma,cl);
  S = D %*% R %*% D
  return(S)
}

linearModelSimulation<- function(n=200,
                                 p=10000,
                                 cl=20,
                                 non_zero_num_cl=5,
                                 non_zero_num_beta=10,
                                 sigma=1,
                                 overlap=1,
                                 rho=0,
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
    beta_x[beta.not.null] <- abs(runif(num_non_zero,0.5,1))*10
    beta.not.null.z <- sample(beta.not.null,size = ceiling(overlap*num_non_zero),replace = FALSE)
    beta_z[beta.not.null.z] <- beta_x[beta.not.null.z]
    y_x <- x %*% beta_x + rnorm(n)
    y_z <- z %*% beta_z + rnorm(n)
    return(list(parameters=parameters,x=x,z=z,beta_x=beta_x,beta_z=beta_z,y_x=y_x,y_z=y_z))
}
sml <- NULL
j <- 1
for(overlap in seq(0,1,0.1)){
  i <- 1
  while(i < 2){
    sml[[j]] <- linearModelSimulation(overlap=overlap,n=100)  
    i <- i + 1
    j <- j + 1
  }
}

load("sml.rda") # very big 15G, takes 5 mins.
res_mat <- NULL
total <- length(sml)
# total <- 5 # takes 1min
for(i in 1:total){
  overlap <- sml[[i]]$parameters$overlap
  x <- sml[[i]]$x
  z <- sml[[i]]$z
  y_x <- sml[[i]]$y_x
  y_z <- sml[[i]]$y_z
  beta.not.null <- which(sml[[i]]$beta_x > 0)
  beta_z <- sml[[i]]$beta_z
  beta_x <- sml[[i]]$beta_x
  mean(abs(cor(x[,beta.not.null],y)))
  mean(abs(cor(x[,which(beta_z > 0)],y)))
  which(sml[[i]]$beta_z > 0)
  method <- 'sis'
  family <- 'gaussian'
  ebic <- FALSE
  ebic.gamma <- 1
  num.select <- dim(x)[2]*0.1
  
  cat(i,"\n") 
  result0 <- screening(x = z,y = y_z,method = 'sis',num.select = num.select,family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  
  result0 <- screening(x = z,y = y_z,method = 'sis',num.select = dim(z)[2],family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  r0 <- sort(result0$screen,decreasing = F,index.return=T)$ix
  #sis
  result1 <- screening(x = x,y = y_x,method = method,num.select = num.select,family = family, ebic = ebic,ebic.gamma = ebic.gamma)
   
  #pool
  result2 <- screening(x = rbind(x,z),y = c(y_x,y_z),method = method,num.select = num.select,family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  
  #ski
  result3 <- SKI(x = x,y = y_x,r0 = r0,method = method,num.select = num.select,family = family,ebic = ebic,ebic.gamma = ebic.gamma)
  
  #ext
  result4 <- screening(x = z,y = y_z,method = method,num.select = num.select,family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  
  sis_tp <- length(which(result1$screen %in% beta.not.null))
  pool_tp <- length(which(result2$screen %in% beta.not.null))
  ski_tp <- length(which(result3$screen %in% beta.not.null))
  ext_tp <- length(which(result4$screen %in% beta.not.null))
  
  alpha_hat <- result3$alpha
  res <- c(sis_tp,pool_tp,ski_tp,ext_tp,alpha_hat,overlap)
  res_mat <- rbind(res_mat,res)
}
res_mat <- cbind(1:dim(res_mat)[1],res_mat)
# draw figure.
rownames(res_mat) <- 1:dim(res_mat)[1]
dat <- as.data.frame(res_mat)
colnames(dat) <- c("iter","sis","pool","ski","ext","alpha","overlap")
library(reshape2)
dat.m <- melt(dat,id.vars='overlap', measure.vars=c("sis","pool","ski","ext"))
dat.m2 <- melt(dat,id.vars='overlap', measure.vars=c("alpha"))
save(dat,file="dat.rda")
save(dat.m,file="dat.m.rda")
save(dat.m2,file="dat.m2.rda")
library(ggplot2)
load("x.rda")
load("dat.rda")
load("dat.m.rda")
load("dat.m2.rda")
p <- ggplot(dat.m) + geom_boxplot(aes(x=as.factor(overlap), y=value, color=variable))
p2 <- ggplot(dat.m2) + geom_boxplot(aes(x=as.factor(overlap), y=value, color=variable))
save(p,file="p.rda")

