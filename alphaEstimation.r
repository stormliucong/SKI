total = 11
for(i in 1:total){
  cat(i,"\n")
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
  family <- 'gaussian'
  ebic <- FALSE
  ebic.gamma <- 1
  num.select <- dim(x)[2]*0.1
  
  result0 <- screening(x = z,y = y_z,method = 'sis',num.select = dim(z)[2],family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  r0 <- sort(result0$screen,decreasing = F,index.return=T)$ix

  result1 <- screening(x = x, y = y_x, method = "sis", num.select=dim(x)[2], family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  r1 <- sort(result$screen,decreasing = F,index.return=T)$ix
  
  
  .alphaEstimation(x=x,y=y_x,r0,r1,num.select,family,method=c("ebic"))
  .alphaEstimation(x=x,y=y_x,r0,r1,num.select,family,method=c("bic"))
  .alphaEstimation(x=x,y=y_x,r0,r1,num.select,family,method=c("deviance"))
}



.alphaEstimation <- function(x,y,r0,r1,num.select.max,family,method=c("ebic","bic","deviance")){
  iter <- 1
  a <- NULL
  t <- NULL
  num_seq <- seq(10,100,by = 10)
  #method=c("ebic")
  for(num in num_seq){
    ix_1 <- which(r0 <= num)
    ix_2 <- which(r1 <= num)
    ix_3 <- sample(1:dim(x)[2],num)
    obj_2 <- cv.glmnet(x[,ix_2],y,family = family)
    fit1 <- glmnet(x[,ix_1],y,family = family,lambda = obj_2$lambda.min)
    fit2 <- glmnet(x[,ix_2],y,family = family,lambda = obj_2$lambda.min)
    fit3 <- glmnet(x[,ix_3],y,family = family,lambda = obj_2$lambda.min)
    if(method == "ebic"){
      v1 <- .ebic(deviance(fit1),num,dim(x)[1],sum(as.numeric(coef(fit1))!=0),1)
      v2 <- .ebic(deviance(fit2),num,dim(x)[1],sum(as.numeric(coef(fit2))!=0),1)
      v3 <- .ebic(deviance(fit3),num,dim(x)[1],sum(as.numeric(coef(fit3))!=0),1)
    }
    if(method == "bic"){
      v1 <- .ebic(deviance(fit1),num,dim(x)[1],sum(as.numeric(coef(fit1))!=0),0)
      v1 <- .ebic(deviance(fit2),num,dim(x)[1],sum(as.numeric(coef(fit2))!=0),0)
      v1 <- .ebic(deviance(fit3),num,dim(x)[1],sum(as.numeric(coef(fit3))!=0),0)
    }
    if(method == "deviance"){
      v1 <- deviance(fit1)
      v2 <- deviance(fit2)
      v3 <- deviance(fit3)
    }
    
    if(v2 >= v3){
      a[iter] <- NA
    }else{
      if(v1 >= v3){
        a[iter] <- 0
      }else{
        a[iter] <- (v1-v3)/(v2-v3)
      }
    }

    not_null <- NULL
    alpha_seq <- seq(0,1,length.out = 6)
    for(alpha in alpha_seq){
      r <- (r0*alpha+r1*(1-alpha))
      rank <- rank(r)
      ix <- which(rank <= num)
      a1=sum(ix %in% beta.not.null)
      not_null <- c(not_null,a1)
    }
    t[iter] <- alpha_seq[which.max(not_null)]
    iter <- iter + 1
  }
  alpha_est <- median(a,na.rm = T)
  true_alpha <- mean(t)
  cat(overlap,"\t",alpha_est,"\t",true_alpha,"\n")
  return(list(overlap=overlap,alpha_est=alpha_est,true_alpha = true_alpha))
}

.ebic <- function(deviance, model.size, sample.size, num.select, ebic.gamma) {
  return (deviance + num.select * (log(sample.size) + 2 * ebic.gamma * log(model.size)))
}
