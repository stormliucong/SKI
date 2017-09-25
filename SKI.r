SKI <- function(x, y, r0, method, num.select, family, ebic, ebic.gamma,cv=FALSE){
  result <- screening(x = x, y = y, method = method, num.select=dim(x)[2], family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  r1 <- sort(result$screen,decreasing = F,index.return=T)$ix
  #current_result <- .alphaEstimation(x = x, y = y, r1 = r1, r0 = r0,alphas = seq(0,1,0.1),num.select.max=num.select.max,family = family)
  #alpha_hat <- current_result$alpha
  alpha_hat <- .alphaEstimation(x = x,y = y,r1 = r1,r0 = r0,num.select.max = 1000,family = family)$alpha
  r <- .combineRank(r0 = r0,r1 = r1,alpha = alpha_hat)
  ix <- which(r <= num.select)
  return(list(alpha=alpha_hat,combined_rank=r,screen=ix))
}

.combineRank <- function(r0,r1,alpha = 0.5){
  for(alpha in seq(0,2,length.out = 10)){
    r <- r0^(alpha/2)*r1^(1-alpha/2)
    rank <- rank(r)
    ix <- which(rank <= num)
    a=sum(ix %in% beta.not.null)
    b=sum(ix_1 %in% beta.not.null)
    c=sum(ix_2 %in% beta.not.null)
    cat("alpha",alpha,"new",a,"ext",b,"int",c,"\n")
  }
  
  return(rank)
}

.alphaEstimation <- function(x,y,r1,r0,num.select.max,family,method=c("ebic","bic","deviance"){
  result <- screening(x = x, y = y, method = method, num.select=dim(x)[2], family = family, ebic = ebic,ebic.gamma = ebic.gamma)
  r1 <- sort(result$screen,decreasing = F,index.return=T)$ix
  iter <- 1
  a <- NULL
  t <- NULL
  for(num in seq(10,num.select.max,length.out = 100)){
    ix_1 <- which(r0 <= num)
    ix_2 <- which(r1 <= num)
    #ix_3 <- sample(1:dim(x)[2],num)
    obj_1 <- cv.glmnet(x[,ix_1],y,family = family)
    obj_2 <- cv.glmnet(x[,ix_2],y,family = family)
    #obj_3 <- cv.glmnet(x[,ix_3],y,family = family)
    fit1 <- glmnet(x[,ix_1],y,family = family,lambda = obj_1$lambda.min)
    fit2 <- glmnet(x[,ix_2],y,family = family,lambda = obj_2$lambda.min)
    #fit3 <- glmnet(x[,ix_3],y,family = family,lambda = obj_3$lambda.1se)
    ebic_1 <- .ebic(deviance(fit1),num,dim(x)[1],sum(as.numeric(coef(fit1))!=0),0)
    ebic_2 <- .ebic(deviance(fit2),num,dim(x)[1],sum(as.numeric(coef(fit2))!=0),)
    ebic_3 <- .ebic(fit2$nulldev,num,dim(x)[1],0,1)
    if(deviance(fit1) >= fit2$nulldev){
      a[iter] <- 0
    }else{
      if(deviance(fit2) >= fit2$nulldev){
        a[iter] <- 1
      }else{
        a[iter] <- (deviance(fit1)-fit2$nulldev)/(deviance(fit2)-fit2$nulldev)
      }
    }
    #a[iter] <- 1-(ebic_2-ebic_1)/(ebic_2-ebic_3)
    
    not_null <- NULL
    alpha_seq <- seq(0,1,length.out = 11)
    for(alpha in alpha_seq){
      r <- r0^alpha*r1^(1-alpha)
      rank <- rank(r)
      ix <- which(rank <= num)
      a1=sum(ix %in% beta.not.null)
      not_null <- c(not_null,a1)
    }
    t[iter] <- alpha_seq[which.max(not_null)]
    iter <- iter + 1
  }
  alpha_est <- mean(a)
  true_alphap <- mean(t)
  
  # for(num in seq(10,num.select.max,length.out = 100)){
  #   ix_1 <- which(r0 <= num)
  #   ix_2 <- which(r1 <= num)
  #   int <- intersect(ix_1,ix_2)
  #   a[iter] <- (length(int)-num^2/dim(x)[2])/num
  #   iter <- iter + 1
  # }
  # 
  # 
  #   
  #   intersect/(num*2-intersect)
  #   cor1 <- mean(abs(cor(x[,ix_1],y)))
  #   cor2 <- mean(abs(cor(x[,ix_2],y)))
  #   cor3 <- abs(cor(x[,ix_random],y))
  #   (cor3 - cor1)/(cor3 - cor2)
  #   obj_1 <- glmnet(x[,ix_1],y,family = family,alpha = 1)
  #   obj_2 <- glmnet(x[,ix_2],y,family = family,alpha = 1)
  #   obj_3 <- glmnet(x[,ix_random],y,family = family,alpha = 1)
  #   ebic_2 <-
  #   ebic_3 <-
  #   bic_1 <- median(deviance(obj_1) + obj_1$df * log(obj_1$nobs))
  #   bic_2 <- median(deviance(obj_2) + obj_2$df * log(obj_2$nobs))
  #   bic_3 <- median(deviance(obj_3) + obj_3$df * log(obj_3$nobs))
  #   (bic_3 - bic_1)/(bic_3 - bic_2)
  #   iter <- iter + 1
  # }
  #   ix_random <- sample(1:dim(x)[2],)
  #   tp_length_1 <- length(which(ix_1 %in% beta.not.null))
  #   tp_length_2 <- length(which(ix_2 %in% beta.not.null))
  #   length(which(ix_random %in% beta.not.null))
  #   obj_1 <- glmnet(x[,ix_1],y,family = family)
  #   obj_2 <- glmnet(x[,ix_2],y,family = family)
  #   obj_3 <- glmnet(x[,])
  #   obj_null <- glmnet(rep(1,dim(x)[1]),y,family = family)
  #   bic_1 <- deviance(obj_1) + obj_1$df * log(obj_1$nobs)
  #   bic_2 <- deviance(obj_2) + obj_2$df * log(obj_2$nobs)
  #   a[i] <- sum(min(bic_1)<min(bic_2))
  #   i <- i + 1
  # }
  # alpha <- sum(a)/100
  return(list(alpha = alpha))
}

# .alphaEstimation <- function(x,y,r1,r0,alphas,num.select.max,family){
#   bic_now <- 1000000
#   num_now <- NULL
#   a_now <- NULL
#   for(a in alphas){
#     for(num in seq(10,num.select.max,length.out = 10)){
#       r_1 <- .combineRank(r0 = r0,r1 = r1,a)
#       ix <- which(r <= num)
#       # Get TP number. Only for test.
#       tp_length <- length(which(ix %in% beta.not.null))
#       obj <- glmnet(x[,ix],y,family = family,alpha = 1)
#       # get loglikehood.
#       d <- deviance(obj)
#       # calculate the BIC.
#       bic <- d + obj$df * log(obj$nobs)
#       cat("a:",a,"num",num,"tp_length:",tp_length,"bic:",min(bic),"\n")
#       if(min(bic) < bic_now){
#         bic_now <- min(bic)
#         num_now <- num
#         a_now <- a
#       }
#     }
#   
#     
#     # dev_ratio <- max(obj$dev.ratio)
#     # dev_ratios <- c(dev_ratios,dev_ratio)
#     # current_alphas <- c(current_alphas,a)
#     # 
#     # if(dev_ratio > current_dev_ratio){
#     #   current_dev_ratio <- dev_ratio
#     #   current_alpha <- a
#     #   cat(a,",",tp_length,",",dev_ratio,"***","\n")
#     # }else{
#     #   cat(a,",",tp_length,",",dev_ratio,"\n")
#     # }
#   }
#   #return(list(alpha_hat = current_alpha, dev_ratio = current_dev_ratio))
#   #return(list(alphas = current_alphas,dev_ratios = dev_ratios))
#   return(list(alpha = a_now, num = num_now, bic = bic_now))
# }

.ebic <- function(deviance, model.size, sample.size, num.select, ebic.gamma) {
  return (deviance + num.select * (log(sample.size) + 2 * ebic.gamma * log(model.size)))
}

