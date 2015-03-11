
library("survival")
library("BB")
library(numDeriv)
library(parallel)
library(splines)
library(survival)
library(orthogonalsplinebasis)
library(MASS)
simwei2 <- function(i, t, l1, l2, l3, b1, b2, b3, a, cov, cen1, cen2){
    expb1 <- exp(sum(cov[i, ] * b1))
    expb2 <- exp(sum(cov[i, ] * b2))
    expb3 <- exp(sum(cov[i, ] * b3))
    lb1 <- l1 * expb1
    lb2 <- l2 * expb2
    lb3 <- l3 * expb3
    p <- lb2/(lb1 + lb2)
    r <- rbinom(1, 1, p)
    u <- runif(3)
    t1 <- ((u[1]^(-t) - 1) / (t * (lb1 + lb2)))^(1/a)
    t2 <- ((u[2]^(-t ) - 1) / (t * (lb1 + lb2)))^(1/a)
    if(r == 0){
        t2 <-  ((u[3]^(- t / (1 + t)) * ( 1 + t * lb1 * t1 ^ a + t * lb2 * t1^a) - ( 1 + t * lb1 * t1 ^ a + t * lb2 * t1^a - t * lb3 * t1 ^ a))/ (t * lb3))^ (1/a)
    }else{
        t1 <- t2 + 3
    }
    c <- runif(1, cen1, cen2) #censoring time
    c(t1, t2, c)
}
simCpRsk <- function(n, p, theta,  lambda1, lambda2, lambda3, kappa, beta1, beta2, beta3, covm = NULL,  cen1, cen2){
    if(is.null(covm)){
        covm <- matrix(rnorm(p * n), n, p) #covariance matrix
    }
    simdata1 <- t(sapply(1 : n, simwei2, theta,  lambda1, lambda2, lambda3, beta1, beta2, beta3, kappa, covm, cen1, cen2))
    d1 <- simdata1[, 1] < simdata1[, 2]&simdata1[, 1] < simdata1[, 3]
    d2 <- simdata1[, 2] < simdata1[, 3]
    y2 <- pmin(simdata1[, 2], simdata1[, 3])
    y1 <- pmin(simdata1[, 1], y2)
    simresp1 <- cbind(y1, d1, y2, d2)
    colnames(simresp1) <- c("y1", "d1", "y2", "d2")
    survData = cbind(simresp1, covm)
    return(survData)
}


    

margpartial4 <- function(paras, resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,   Av1, Av2 ){
    
    pl1 <- ncol(A1)
    pl2 <- ncol(A2)
    pl3 <- ncol(A32)
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    sp1 <- paras[ (3 * p + 1): (3 * p + pl1)]
    sp2 <- paras[ (3 * p + pl1 + 1): (3 * p + pl1 + pl2)]
    sp3 <- paras[ (3 * p + pl1 + pl2 + 1): (3 * p + pl1 + pl2 + pl3)]
    theta <- paras[  (3 * p + pl1 + pl2  +pl3) + 1]
#   mtheta <- c(mtheta, theta)
 #   mtheta <<- mtheta[max(length(mtheta)- 500, 1) : length(mtheta)]
    Lambda1 <- A1 %*% sp1
    Lambda2 <- A2  %*% sp2
    Lambda32 <- A32 %*% sp3
    Lambda31 <- A31 %*% sp3
    ## lambda1 <- pmax(predict(smooth.spline(resp[, "y1"], Lambda1), resp[ind1, "y1"], deriv= 1)$y, 0)
    ## lambda2 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda2), resp[ind2, "y2"], deriv= 1)$y, 0)
    ## lambda3 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda32), resp[ind3, "y2"], deriv= 1)$y, 0)
    ## if(sum(lambda1<0) > 0|sum(lambda2<0) > 0 |sum(lambda3<0) > 0){
    ##     browser()
    ## }
    lambda1 <- dA1 %*% sp1
    lambda2 <- dA2 %*% sp2
    lambda3 <- dA3 %*% sp3
    
  
    A <- (Lambda1 -  Av1 %*% sp1)  * exp(cov %*% beta1) +  (Lambda2 -  Av2 %*% sp2)  * exp(cov %*% beta2) +  (Lambda32 - Lambda31)  * exp(cov %*% beta3)
    vl <- c( (lambda1 ), (lambda2), (lambda3) )#c(lambda1, lambda2, lambda3)#c(lambda1, lambda2, lambda3)##
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16)))) 
    if(sum(is.nan(log(vl + 1e-16))) > 0){
        browser()
    }else{
        return(-res)
    }
    

    
}


margpartial4sum <- function(paras,  resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,   Av1, Av2 ){
    
    pl1 <- ncol(A1)
    pl2 <- ncol(A2)
    pl3 <- ncol(A32)
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    sp1 <- paras[ (3 * p + 1): (3 * p + pl1)]
    sp2 <- paras[ (3 * p + pl1 + 1): (3 * p + pl1 + pl2)]
    sp3 <- paras[ (3 * p + pl1 + pl2 + 1): (3 * p + pl1 + pl2 + pl3)]
    theta <- paras[  (3 * p + pl1 + pl2 + pl3) + 1]
#   mtheta <- c(mtheta, theta)
 #   mtheta <<- mtheta[max(length(mtheta)- 500, 1) : length(mtheta)]
    Lambda1 <- A1 %*% sp1
    Lambda2 <- A2  %*% sp2
    Lambda32 <- A32 %*% sp3
    Lambda31 <- A31 %*% sp3
    ## lambda1 <- pmax(predict(smooth.spline(resp[, "y1"], Lambda1), resp[ind1, "y1"], deriv= 1)$y, 0)
    ## lambda2 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda2), resp[ind2, "y2"], deriv= 1)$y, 0)
    ## lambda3 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda32), resp[ind3, "y2"], deriv= 1)$y, 0)
    ## if(sum(lambda1<0) > 0|sum(lambda2<0) > 0 |sum(lambda3<0) > 0){
    ##     browser()
    ## }
    lambda1 <- dA1 %*% sp1
    lambda2 <- dA2 %*% sp2
    lambda3 <- dA3 %*% sp3
    
  
    A <- (Lambda1 -  Av1 %*% sp1)  * exp(cov %*% beta1) +  (Lambda2 -  Av2 %*% sp2)  * exp(cov %*% beta2) +  (Lambda32 - Lambda31)  * exp(cov %*% beta3)
    vl <- c( (lambda1 ), (lambda2), (lambda3) )#c(lambda1, lambda2, lambda3)#c(lambda1, lambda2, lambda3)##
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16)))) 
    if(sum(is.nan(log(vl + 1e-16))) > 0){
        browser()
    }else{
        return(-res)
    }
    

    
}



margpartial4nosum <- function(paras,  resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,   Av1, Av2 ){
    
    pl1 <- ncol(A1)
    pl2 <- ncol(A2)
    pl3 <- ncol(A32)
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    sp1 <- paras[ (3 * p + 1): (3 * p + pl1)]
    sp2 <- paras[ (3 * p + pl1 + 1): (3 * p + pl1 + pl2)]
    sp3 <- paras[ (3 * p + pl1 + pl2 + 1): (3 * p + pl1 + pl2 + pl3)]
    theta <- paras[  (3 * p + pl1 + pl2 + pl3) + 1]
#   mtheta <- c(mtheta, theta)
 #   mtheta <<- mtheta[max(length(mtheta)- 500, 1) : length(mtheta)]
    Lambda1 <- A1 %*% sp1
    Lambda2 <- A2  %*% sp2
    Lambda32 <- A32 %*% sp3
    Lambda31 <- A31 %*% sp3
    ## lambda1 <- pmax(predict(smooth.spline(resp[, "y1"], Lambda1), resp[ind1, "y1"], deriv= 1)$y, 0)
    ## lambda2 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda2), resp[ind2, "y2"], deriv= 1)$y, 0)
    ## lambda3 <- pmax(predict(smooth.spline(resp[, "y2"], Lambda32), resp[ind3, "y2"], deriv= 1)$y, 0)
    ## if(sum(lambda1<0) > 0|sum(lambda2<0) > 0 |sum(lambda3<0) > 0){
    ##     browser()
    ## }
    lambda1 <- dA1 %*% sp1
    lambda2 <- dA2 %*% sp2
    lambda3 <- dA3 %*% sp3
    mvl1 <- mvl2 <- mvl3 <- rep(1, n)
    mcov1 <- mcov2 <- mcov3 <- matrix(0, n, p)
    mvl1[ind1] <- lambda1 + 1e-16
    mvl2[ind2] <- lambda2 + 1e-16
    mvl3[ind3] <- lambda3 + 1e-16
    mcov1[ind1, ] <- cov1
    mcov2[ind2, ] <- cov2
    mcov3[ind3, ] <- cov3
  
    A <- (Lambda1 -  Av1 %*% sp1)  * exp(cov %*% beta1) +  (Lambda2 -  Av2 %*% sp2)  * exp(cov %*% beta2) +  (Lambda32 - Lambda31)  * exp(cov %*% beta3)
    vl <- c( (lambda1 ), (lambda2), (lambda3) )#c(lambda1, lambda2, lambda3)#c(lambda1, lambda2, lambda3)##
    sumbb <- (matrix(mcov1, ncol = p)%*%  matrix(beta1))  + (matrix(mcov2, ncol = p)%*%  matrix(beta2))+ (matrix(mcov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try(((resp[, 1] * resp[, 2]) * log(theta + 1)  -(B * log(1 + theta* A))+  sumbb + log(mvl1) + log(mvl2) + log(mvl3))) 
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-res)
    }
    

    
}






iniestreal4 <- function(theta, bb, upper, lower,  resp, cov, nk, ord, var, ctrl = list(M = 10, maxit = 30000, maxfeval = 1e6, trace = FALSE)){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")

     nml <- max(resp[, 3:5])
     resp[, 3:5] <- resp[, 3:5]
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
      knots1 <- quantile(resp[, c(3:5)], seq(0, 1, length.out = nk))#seq(quantile(resp[ind1, c(3)], 0.2), quantile(resp[ind1, c(3)], 0.8), length.out = nk)
knots1 <- expand.knots(knots1, ord)

     knots2 <-quantile(resp[, c(3:5)], seq(0, 1, length.out = nk)) #seq(quantile(resp[ind2, c(3)], 0.2), quantile(resp[ind2, c(3)], 0.8), length.out = nk)#quantile(resp[ind2, c(3,  5 )], seq(0.2, 0.8, length = nk))  #
knots2 <- expand.knots(knots2, ord)
     
     knots3 <- quantile(resp[, c(3:5)], seq(0, 1, length.out = nk))
knots3 <- expand.knots(knots3, ord)
     uniqtime <- unique(as.vector(resp[, 3:5]))
     uniqtime <- uniqtime[order(uniqtime)]
     Bs1 <- SplineBasis(knots1, ord, TRUE)
     Bs2 <- SplineBasis(knots2, ord, TRUE)
     Bs3 <- SplineBasis(knots3, ord, TRUE)
     dA1 <- evaluate(Bs1, resp[ind1, "y1"])
     dA2 <- evaluate(Bs2, resp[ind2, "y1"])
     dA3 <- evaluate(Bs3, resp[ind3, "y2"])
     
     mA1 <- evaluate(integrate(Bs1),uniqtime)
     mA2 <- evaluate(integrate(Bs2),uniqtime)
     mA3 <- evaluate(integrate(Bs3),uniqtime)
     rownames(mA1) <- as.character(uniqtime)
     rownames(mA2) <- as.character(uniqtime)
     rownames(mA3) <- as.character(uniqtime)
    

     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
   
     
     subgix1 <- (d1 == 0)
     respsub1 <- resp[subgix1, ]
     respsub2 <- resp[ind1, ]
     surv1 <- (coxph(Surv(resp[, "y1"], resp[, "d1"]) ~ ., as.data.frame(cov)))
     surv2 <- (coxph(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~., as.data.frame(cov[subgix1, ] )))
     surv3 <- (coxph(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~., as.data.frame(cov[ind1, ])))
     bz10 <- basehaz(surv1)
     bz20 <- basehaz(surv2)
     bz30 <- basehaz(surv3)
     cov1 <- cov[ind1, ]
     cov2 <- cov[ind2, ]
     cov3 <- cov[ind3, ]
     A1 <- mA1[as.character(bz10$time), ]
     A2 <- mA2[as.character(bz20$time), ]
     A3 <- mA3[as.character(bz30$time), ]
     Av1 <- mA1[as.character(resp[, "v"]), ]
     Av2 <- mA2[as.character(resp[, "v"]), ]
     
     sp1 <- pmax(ginv(t(A1) %*% (A1)) %*% t(A1) %*% bz10$hazard, 0)
     sp2 <- pmax(ginv(t(A2) %*% (A2)) %*% t(A2) %*% bz20$hazard, 0)
     sp3 <- pmax(ginv(t(A3) %*% (A3)) %*% t(A3) %*% bz30$hazard, 0)
     A1 <- mA1[as.character(resp[, "y1"]), ]
     A2 <- mA2[as.character(resp[, "y1"]), ]
     A32 <- mA3[as.character(resp[, "y2"]), ]
     A31 <- mA3[as.character(resp[, "y1"]), ]
     pl1 <- ncol(mA1)
     pl2 <- ncol(mA2)
     pl3 <- ncol(mA3)
 #   mtheta <<- 0
     
    
     if(var == TRUE){
         #c(bb[1 : (3 * p)], bb[length(bb)]), bb[(3 * p + 1) : (3 * p + pl1)], bb[(3 * p + pl1 + 1) : (3 * p + pl1 + pl2)], bb[(3 * p +pl1 + pl2 +  1) : (3 * p + pl1 + pl2 + pl3)],
         #c(bb[1 : (3 * p)], bb[length(bb)]), bb[(3 * p + 1) : (3 * p + pl1)], bb[(3 * p + pl1 + 1) : (3 * p + pl1 + pl2)], bb[(3 * p +pl1 + pl2 +  1) : (3 * p + pl1 + pl2 + pl3)],
         #browser()
     D1 <- jacobian(margpartial4nosum, bb,   method="Richardson", method.args=list(), resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,  Av1, Av2)
     D2 <- numDeriv::hessian(margpartial4, bb,  method="Richardson", method.args=list(), resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,  Av1, Av2) 
     gD2 <- ginv(D2)
     var <- t(gD2) %*% t(D1) %*% D1 %*% (gD2)
     res <- NULL
     return(var)
 }else{
      res <- spg(c(beta1, beta2, beta3, sp1, sp2, sp3, theta), margpartial4, gr = NULL, method=2,  lower = lower,  upper = upper,  project=NULL, projectArgs=NULL, control=ctrl,  quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,  Av1, Av2)
     var <- NULL
      return(list(res, Bs1, Bs2, Bs3, p = p))
 }
   
 }



resfun <- function(survData,  iniv, lower, upper, nk, ord,scale, var,  ctrl = list(M = 10, maxit = 30000, maxfeval = 1e6, trace = FALSE)){
    if(var != 1){
        lpara <- length(iniv)
        theta <- iniv[lpara]
        bb <- iniv[1 : (lpara - 1)]
    }else{
            bb <- iniv
        }
    survData[, c(1, 3)] <- survData[, c(1, 3)]/scale#/ max(survData[, c(1, 3)])
    np <- ncol(survData)
    n <- nrow(survData)
    y1 <- pmin(survData[, 1], survData[, 3])
    y2 <- survData[, 3]
    d1 <- survData[, 2]
    d2 <- survData[, 4]

    v <- rep(0, n)
    covmy <- matrix(survData[, (5 : np)], ncol =  np - 4)
   

       
resp <- cbind(d1, d2, y1, y2, v)
colnames(resp) <- cbind("d1", "d2", "y1", "y2", "v")
res <-  iniestreal4(theta, bb, upper, lower,  resp, covmy, nk, ord, var, ctrl) #apply(mres, 2, median)
 
}

FrqSp <- function(survData, startValues, lower, upper, nknots, ord, scale,  ctrl = list(M = 10, maxit = 30000, maxfeval = 1e6, trace = FALSE)){
    tres <- resfun(survData, startValues, lower, upper, nknots, ord, scale, 0,  ctrl)
    tvar <- resfun(survData, tres[[1]][[1]], lower, upper, nknots, ord, scale, 1)
    res <- list(tres, tvar, scale)
    class(res) <- "FrqSp"
    return(res)
}

summary.FrqSp <- function(objects){
    p<- objects[[1]]$p
    beta1 <- objects[[1]][[1]]$par[ 1:(p)]
    beta2 <- objects[[1]][[1]]$par[ (p+ 1):(2 * p)]
    beta3 <- objects[[1]][[1]]$par[ (2 * p+ 1):(3 * p)]
    theta <- objects[[1]][[1]]$par[length(objects[[1]][[1]]$par)]
    msd <- sqrt(diag(objects[[2]]))
    sbeta1 <- msd[ 1:(p)]
    sbeta2 <- msd[ (p + 1 ): (2 * p)]
    sbeta3 <- msd[(2 * p + 1) : (3 *p)]
    rbeta1 <- cbind(beta1, sbeta1)
    rbeta2 <- cbind(beta2, sbeta2)
    rbeta3 <- cbind(beta3, sbeta3)
    rtheta <- cbind(theta, msd[length(objects[[1]][[1]]$par)])
    colnames(rbeta1) <- colnames(rbeta2) <- colnames(rbeta3) <- c("beta", "sd")
    colnames(rtheta) <- c("theta", "sd")
    list(beta1 = rbeta1, beta2 = rbeta2, beta3 = rbeta3, theta = rtheta)
}
plot.FrqSp <- function(objects, lty = c(1, 2, 3), ylim = NULL){
    p <- objects[[1]]$p
    Bs1 <- objects[[1]][[2]]
    Bs2 <- objects[[1]][[3]]
    Bs3 <- objects[[1]][[4]]
   
    tm <- (seq(0, 1, 0.01) * objects[[3]])
    dA1 <- evaluate(integrate(Bs1), seq(0, 1, 0.01))
    dA2 <- evaluate(integrate(Bs2), seq(0, 1, 0.01))
    dA3 <- evaluate(integrate(Bs3), seq(0, 1, 0.01))
    lp <- ncol(dA1) 
    sp1 <- objects[[1]][[1]]$par[(3 * p + 1):( 3 * p + 1 + lp - 1)]
    sp2 <- objects[[1]][[1]]$par[( 3 * p + 1 + lp ):( 3 * p + 1 + 2 * lp - 1)]
    sp3 <- objects[[1]][[1]]$par[( 3 * p + 1 + 2 * lp ):( 3 * p + 1 + 3 * lp - 1)]
    bspline1 <- ( (dA1 %*% sp1) ) 
    bspline2 <- ( (dA2 %*% sp2 ))
    bspline3 <- ( (dA3 %*% sp3))
    plot(bspline1 ~ tm, type = "l", xlab = "Time to event", ylab = "Cumulative hazard", ylim = ylim, lty = lty[1])
    lines(bspline2 ~ tm, type = "l", ylab = "Cumulative hazard", ylim = ylim, lty = lty[2] )
    lines(bspline3 ~ tm, type = "l",  ylab = "Cumulative hazard", ylim = ylim, lty = lty[3])
   legend("topleft", c("nonterminal event", "terminal event w/o nonterminal event", "terminal event with nonterminal event"), lty = lty, cex = 0.8)
        
    
    
}
        
