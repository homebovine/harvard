
library("survival")
library("BB")
library(numDeriv)
library(rootSolve)
library(parallel)
library(splines)
library(survival)
library(orthogonalsplinebasis)
library(MASS)
library(ICsurv)


    
margpartial <- function(paras, resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,   Av1, Av2 ){
    
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
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16)))/n)
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-res)
    }
    

    
}






iniestreal <- function(theta, bb,  resp, cov){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")
   #  resp[, 3:5] <- exp(resp[, 3:5])/(1 + exp(resp[, 3:5]))
     nml <- max(resp[, 3:5])
     resp[, 3:5] <- resp[, 3:5]
     ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
      knots1 <- quantile(resp[ind1, c(3)], seq(0, 1, length.out = nk))#seq(quantile(resp[ind1, c(3)], 0.2), quantile(resp[ind1, c(3)], 0.8), length.out = nk)
knots1 <- c(unique(max(min(knots1) - 3 * mean(diff(knots1)), 0), max(min(knots1) - 2* mean(diff(knots1)), 0), max(min(knots1) - mean(diff(knots1)), 0)), knots1, max(knots1) + mean(diff(knots1)), max(knots1) + 2 * mean(diff(knots1)), max(knots1) + 3 * mean(diff(knots1)))

     knots2 <-quantile(resp[ind2, c(3)], seq(0, 1, length.out = nk)) #seq(quantile(resp[ind2, c(3)], 0.2), quantile(resp[ind2, c(3)], 0.8), length.out = nk)#quantile(resp[ind2, c(3,  5 )], seq(0.2, 0.8, length = nk))  #
knots2 <- c(unique(max(min(knots2) - 3 * mean(diff(knots2)), 0), max(min(knots2) - 2* mean(diff(knots2)), 0), max(min(knots2) - mean(diff(knots2)), 0)), knots2, max(knots2) + mean(diff(knots2)), max(knots2) + 2 * mean(diff(knots2)), max(knots2) + 3 * mean(diff(knots2)))
     
     knots3 <- quantile(resp[ind3, c(4)], seq(0, 1, length.out = nk))
knots3 <- c(unique(max(min(knots3) - 3 * mean(diff(knots3)), 0), max(min(knots3) - 2* mean(diff(knots3)), 0), max(min(knots3) - mean(diff(knots3)), 0)), knots3, max(knots3) + mean(diff(knots3)), max(knots3) + 2 * mean(diff(knots3)), max(knots3) + 3 * mean(diff(knots3)))#
     uniqtime <- unique(as.vector(resp[, 3:5]))
     uniqtime <- uniqtime[order(uniqtime)]
#Bs1 <- SplineBasis(knots1, ord, TRUE)
  #    Bs2 <- SplineBasis(knots2, ord, TRUE)
#      Bs3 <- SplineBasis(knots3, ord, TRUE)

     ## Bs1 <- SplineBasis(knots1, ord, TRUE)
##      Bs2 <- SplineBasis(knots2, ord, TRUE)
##      Bs3 <- SplineBasis(knots3, ord, TRUE)
##      mdA1 <- evaluate(Bs1, uniqtime)
##      #mdA2 <- evaluate(Bs2, uniqtime)
##      #mdA3 <- evaluate(Bs3, uniqtime)
## #     mdA1[is.na(mdA1)] <- 0
##      mdA3 <- mdA2 <- mdA1
     rnA <- as.character(uniqtime)
     
  #   rownames(dA) <- rnA
     lt <- length(uniqtime)
     mA1 <- t(Ispline(uniqtime, ord, knots1))
     mA2 <- t(Ispline(uniqtime, ord, knots2))
     mA3 <- t(Ispline(uniqtime, ord, knots3))
     
     rownames(mA1) <- rnA
     rownames(mA2) <- rnA
     rownames(mA3) <- rnA
     ## rownames(mdA1) <- rnA
     ## rownames(mdA2) <- rnA
     ## rownames(mdA3) <- rnA
     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     
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
     ## dA1 <- dA[as.character(resp[, "y1"]), ]
     ## dA2 <- dA[as.character(resp[, "y2"]), ]
     ## dA3 <- dA2
     sp1 <- pmax(ginv(t(A1) %*% (A1)) %*% t(A1) %*% bz10$hazard, 0.01)
     sp2 <- pmax(ginv(t(A2) %*% (A2)) %*% t(A2) %*% bz20$hazard, 0.01)
     sp3 <- pmax(ginv(t(A3) %*% (A3)) %*% t(A3) %*% bz30$hazard, 0.01)
     A1 <- mA1[as.character(resp[, "y1"]), ]
     A2 <- mA2[as.character(resp[, "y1"]), ]
     A32 <- mA3[as.character(resp[, "y2"]), ]
     A31 <- mA3[as.character(resp[, "y1"]), ]
     pl1 <- ncol(mA1)
     pl2 <- ncol(mA2)
     pl3 <- ncol(mA3)
     o <- order(resp[ind1, "y1"])
     dA1 <- apply(rbind(rep(0,pl1) ,mA1[as.character(resp[ind1, "y1"][o]), ]), 2, diff)
   
     
     o <- order(resp[ind2, "y1"])
     dA2 <- apply(rbind(rep(0,pl2) ,mA2[as.character(resp[ind2, "y1"][o]), ]), 2, diff)

     
     o <- order(resp[ind3, "y2"])
     dA3 <- apply(rbind(rep(0,pl3) ,mA3[as.character(resp[ind3, "y2"][o]), ]), 2, diff)
     
     
     res <- spg(c(beta1, beta2, beta3, sp1, sp2, sp3, theta), margpartial, gr = NULL, method=2,  lower = c(rep(-1, 3 * p), rep(0, pl1 + pl2 + pl3), 0.01), upper = c(rep(10, 3 * p), rep(10, pl1 + pl2 + pl3), 1), project=NULL, projectArgs=NULL, control=list(M = 10, maxit = 5000), quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,  Av1, Av2)
    res
 }


#rescen05084 <- FrqID( survData, rep(0, 9), stheta = c(0.5, 0.84), tol = 1e-6,  ltr = T, step = 0.02,ncores = 10,   verbose =2)







    
margpartial1 <- function(paras, resp, cov, n, p, cov1, cov2, cov3, ind1, ind2, ind3, o1, o2,  knots1, knots2, knots3, uniqtime1, rname1, dtime1, uniqtime2, rname2, dtime2){
   
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    dA1 <- paras[ (3 * p + 1): (3 * p + pl)]
    dA2 <- paras[ (3 * p + pl + 1): (3 * p + 2 * pl)]
    dA3 <- paras[ (3 * p + 2 * pl + 1): (3 * p + 3 * pl)]
    theta <-  paras[(3 * p + 3 * pl) + 1]
#    f <- paras[(3 * p + 3 * pl) + 2]
    dA1 <- (dA1)
    dA2 <- (dA2)
    dA3 <- (dA3) 
    lambda1 <- (as.matrix(approxfun(knots1, dA1, method = "constant", f= 1/2, rule = 2)(uniqtime1)))
                         
    lambda2 <- (as.matrix(approxfun(knots2, dA2, method = "constant", f= 1/2, rule = 2)(uniqtime1)))
    lambda3 <- (as.matrix(approxfun(knots3, dA3, method = "constant", f= 1/2, rule = 2)(uniqtime2)))
    if(sum(is.na(c(lambda1, lambda2, lambda3))) >0){
        browser()
    }
    rownames(lambda1) <- rname1
    rownames(lambda2) <- rname1 
    rownames(lambda3) <- rname2 
   
    Lambda1 <- as.matrix(cumsum(lambda1 * dtime1 ))
    Lambda2 <- as.matrix(cumsum(lambda2 * dtime1 ))
    Lambda3 <- as.matrix(cumsum(lambda3 * dtime2 ))
    rownames(Lambda1) <- rname1
    rownames(Lambda2) <- rname1
    rownames(Lambda3) <- rname2
    
    A <- (Lambda1[ as.character(resp[, "y1"]), ] - Lambda1[as.character(resp[, "v"]), ] )  * exp(cov %*% beta1) + (Lambda2[ as.character(resp[, "y1"]), ] - Lambda2[as.character(resp[, "v"]), ] )  * exp(cov %*% beta2) +  (Lambda3[ as.character(resp[, "y2"]), ] - Lambda3[ as.character(resp[, "y1"]), ])  * exp(cov %*% beta3)
    vl <-c( ((lambda1* dtime1)[as.character(resp[ind1, "y1"]), ] ), ((lambda2 * dtime1)[as.character(resp[ind2, "y1"]), ]), ((lambda3 * dtime2)[as.character(resp[ind3, "y2"]), ]) )#c(lambda1, lambda2, lambda3)# c(lambda1, lambda2, lambda3)##
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16))))
    #print(theta)
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-res)
    }
    

    
}






iniestreal1 <- function(theta, bb,  resp, cov){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")
     nml <- max(resp[, 3:5])
     #resp[, 3:5] <- exp(resp[, 3:5])/(1 + exp(resp[, 3:5)])
     o1 <- order(resp[, "y1"])
     o2 <- order(resp[, "y2"])
     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
     knots1 <- seq(quantile(resp[ind1, c(3)], 0.2), quantile(resp[ind1, c(3)], 0.8), length.out = nk)
     knots2 <- seq(quantile(resp[ind2, c(4)], 0.2), quantile(resp[ind2, c(4)], 0.8), length.out = nk)#quantile(resp[ind2, c(3,  5 )], seq(0.2, 0.8, length = nk))  #
     knots3 <- seq(quantile(resp[ind3, c(4)],0.2 ), quantile(resp[ind3, c(4)], 0.8), length.out = nk)#quantile(resp[ind3, c(3,  4)], seq(0.2, 0.8, length = nk)) #quantile(resp[, c(3:5)], seq(0, 1, length = nk))#
     uniqtime1 <- unique(as.vector(resp[, c(3, 5)]))
     uniqtime1 <- uniqtime1[order(uniqtime1)]
     rname1 <- as.character(uniqtime1)
     dtime1 <- as.matrix(c(1, diff(uniqtime1)))
     rownames(dtime1) <- rname1
     uniqtime2 <- unique(as.vector(resp[, c(3, 4)]))
     uniqtime2 <- uniqtime2[order(uniqtime2)]
     rname2 <- as.character(uniqtime2)
     dtime2 <- as.matrix(c(1, diff(uniqtime2)))
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
     hbz10 <- bz10
     hbz10$hazard <- c(bz10$hazard[1], diff(bz10$hazard))/c(1, diff(bz10$time))
     hbz20 <- bz20
     hbz20$hazard <- c(bz20$hazard[1], diff(bz20$hazard))/c(1, diff(bz20$time))
     hbz30 <- bz30
     hbz30$hazard <- c(bz30$hazard[1], diff(bz30$hazard))/c(1, diff(bz30$time))
     dA1 <- (approxfun(hbz10$time[hbz10$hazard > 0], (hbz10$hazard[hbz10$hazard > 0]), method = "constant", f= 1/2, rule = 2)(knots1))
     dA2 <- (approxfun(hbz20$time[hbz20$hazard > 0], (hbz20$hazard[hbz20$hazard > 0]), method = "constant", f= 1/2, rule = 2)(knots2))
     dA3 <- (approxfun(hbz30$time[hbz30$hazard > 0], (hbz30$hazard[hbz30$hazard > 0]), method = "constant", f= 1/2, rule = 2)(knots3))
  
     res <- spg(c(beta1, beta2, beta3, dA1, dA2, dA3, theta), margpartial1, gr = NULL, method=2,  lower = c(rep(-1, 3 * p), rep(0.01, 3 * pl), 0.01), upper = c(rep(10, 3 * p), rep(10, 3 * pl), 10), project=NULL, projectArgs=NULL, control=list(M = 10, ftol = 1.e-10, maxit = 1000), quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3,  ind1, ind2, ind3, o1, o2, knots1, knots2, knots3,  uniqtime1, rname1, dtime1, uniqtime2, rname2, dtime2)
    res
 }



margpartial2 <- function(paras, resp, cov, n, p, cov1, cov2, cov3, ind1, ind2, ind3, o1, o2,  knots1, knots2, knots3, pl1, pl2, pl3, uniqtime1, rname1, dtime1, intgtime1, intgname1){
   
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    dA1 <- paras[ (3 * p + 1): (3 * p + pl1- 1)]
    dA2 <- paras[ (3 * p + pl1- 1 + 1): (3 * p + pl1- 1 + pl2- 1)]
    dA3 <- paras[ (3 * p + pl1- 1 + pl2- 1 + 1): (3 * p + pl1 - 1+ pl2- 1 + pl3- 1)]
    theta <-  paras[(3 * p + pl1 + pl2 + pl3- 3) + 1]
    #w <- 0#paras[((3 * p + pl1 + pl2 + pl3) + 2)]
    
   # cumA1 <- cumsum(dA1)
    #cumA2 <- cumsum(dA2)
    #cumA3 <- cumsum(dA3)
  #  f <- paras[(3 * p + 3 * pl) + 2]
    dA1 <- c(0, dA1)
    dA2 <- c(0, dA2)
    dA3 <- c(0, dA3)
   # if(as.integer(dA1[1]) != 1){browser()}
    ## lambda11 <- as.matrix(stepfun(x = knots1, y = dA1, f= 1)(uniqtime1))
    ## lambda21 <- as.matrix(stepfun(x = knots2, y = dA2, f= 1)(uniqtime1))
    ## lambda31 <- as.matrix(stepfun(x = knots3, y = dA3, f= 1)(uniqtime2))

    # lambda12 <- as.matrix(stepfun(x = knots1, y = dA1, f= w)(intgtime1))
    # lambda22 <- as.matrix(stepfun(x = knots2, y = dA2, f= w)(intgtime1))
    # lambda32 <- as.matrix(stepfun(x = knots3, y = dA3, f= w)(intgtime1))

    ## lambda13 <- as.matrix(stepfun(x = knots1, y = dA1, f= 3/4)(uniqtime1))
    ## lambda23 <- as.matrix(stepfun(x = knots2, y = dA2, f= 3/4)(uniqtime1))
    ## lambda33 <- as.matrix(stepfun(x = knots3, y = dA3, f= 3/4)(uniqtime2))

    ## lambda14 <- as.matrix(stepfun(x = knots1, y = dA1, f= 1/4)(uniqtime1))
    ## lambda24 <- as.matrix(stepfun(x = knots2, y = dA2, f= 1/4)(uniqtime1))
    ## lambda34 <- as.matrix(stepfun(x = knots3, y = dA3, f= 1/4)(uniqtime2))

    lambda1 <- as.matrix(stepfun(x = knots1, y = dA1, f= 0)(intgtime1))
    lambda2 <- as.matrix(stepfun(x = knots2, y = dA2, f= 0)(intgtime1))
    lambda3 <- as.matrix(stepfun(x = knots3, y = dA3, f= 0)(intgtime1))
    if(sum(is.na(c(lambda1, lambda2, lambda3))) >0){
        browser()
    }
    rownames(lambda1) <- intgname1
    rownames(lambda2) <- intgname1 
    rownames(lambda3) <- intgname1 
    intgral <- function(y, knots, dA){
        stats::integrate(stepfun(x = knots, y = dA, f= w), 0, y)$value
    }
   
    Lambda1 <- as.matrix(cumsum(lambda1 * dtime1 )) 
    Lambda2 <- as.matrix(cumsum(lambda2 * dtime1 )) 
    Lambda3 <- as.matrix(cumsum(lambda3 * dtime1 )) 
    rownames(Lambda1) <- intgname1
    rownames(Lambda2) <- intgname1
    rownames(Lambda3) <- intgname1
    
    A <- (Lambda1[ as.character(resp[, "y1"]), ] - Lambda1[as.character(resp[, "v"]), ] )  * exp(cov %*% beta1) + (Lambda2[ as.character(resp[, "y1"]), ] - Lambda2[as.character(resp[, "v"]), ] )  * exp(cov %*% beta2) +  (Lambda3[ as.character(resp[, "y2"]), ] - Lambda3[ as.character(resp[, "y1"]), ])  * exp(cov %*% beta3)
    vl <-c( ((lambda1)[as.character(resp[ind1, "y1"]), ] ), ((lambda2  )[as.character(resp[ind2, "y1"]), ]), ((lambda3  )[as.character(resp[ind3, "y2"]), ]) )#c(lambda1, lambda2, lambda3)# c(lambda1, lambda2, lambda3)##
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16)))/n)
    #print(theta)
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-
    
}






iniestreal2 <- function(theta, bb,  resp, cov, knots1=NULL, knots2 = NULL, knots3 = NULL){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")
     nml <- max(resp[, 3:5])
    # resp[, 3:5] <- exp(resp[, 3:5])/(1 + exp(resp[, 3:5]))
     o1 <- order(resp[, "y1"])
     o2 <- order(resp[, "y2"])
     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
     cov1 <- cov[ind1, ]
     cov2 <- cov[ind2, ]
     cov3 <- cov[ind3, ]
     if (is.null(knots1))
     knots1 <- quantile(resp[ind1, c(3)], seq(0, 1, length.out = nk))#seq(quantile(resp[ind1, c(3)], 0.2), quantile(resp[ind1, c(3)], 0.8), length.out = nk)
     if (is.null(knots2))
     knots2 <-quantile(resp[ind2, c(3)], seq(0, 1, length.out = nk)) #seq(quantile(resp[ind2, c(3)], 0.2), quantile(resp[ind2, c(3)], 0.8), length.out = nk)#quantile(resp[ind2, c(3,  5 )], seq(0.2, 0.8, length = nk))  #
     if (is.null(knots3))
     knots3 <- quantile(resp[ind3, c(4)], seq(0, 1, length.out = nk))#seq(quantile(resp[ind3, c(4)],0.2 ), quantile(resp[ind3, c(4)], 0.8), length.out = nk)#quantile(resp[ind3, c(3,  4)], seq(0.2, 0.8, length = nk)) #quantile(resp[, c(3:5)], seq(0, 1, length = nk))#
     uniqtime1 <- unique(as.vector(resp[, c(3, 4, 5)]))
     uniqtime1 <- uniqtime1[order(uniqtime1)]
     intervals <- cbind(uniqtime1, c(uniqtime1[-1], 100))[1:(length(uniqtime1)-1), ]
     fill <- function(interval){
         runif(5, interval[1],interval[2])
     }
     intgtime1 <- c(uniqtime1, apply(intervals, 1, fill))
     intgtime1 <- unique(intgtime1)
     intgtime1 <- intgtime1[order(intgtime1)]
     rname1 <- as.character(uniqtime1)
     intgname1 <- as.character(intgtime1)
     dtime1 <- as.matrix(c(0, diff(intgtime1)))
     ## uniqtime2 <- unique(as.vector(resp[, c(3, 4, 5)]))
     ## uniqtime2 <- uniqtime2[order(uniqtime2)]
     ## rname2 <- as.character(uniqtime2)
     ## dtime2 <- as.matrix(c(1, diff(uniqtime2)))
     subgix1 <- (d1 == 0)
     respsub1 <- resp[subgix1, ]
     respsub2 <- resp[ind1, ]
    
    
     pl1 <- length(knots1) + 1
     pl2 <- length(knots2) + 1
     pl3 <- length(knots3) + 1
     dA1 <-  rep(1, pl1- 1 )
     dA2 <- rep(1, pl1- 1 )
     dA3 <- rep(1, pl1- 1 )
     w <- c(1/2)
     res <- spg(c(beta1, beta2, beta3, dA1, dA2, dA3, theta), margpartial2, gr = NULL, method=2,  lower = c(rep(-1, 3 * p), rep(0, pl1 + pl2+ pl3 - 3), 0.01), upper = c(rep(100, 3 * p), rep(100,  pl1 + pl2+ pl3- 3), 10), project=NULL, projectArgs=NULL, control=list(M = 10, ftol = 1.e-10, maxit = 4000, gtol= 1e-5), quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3,  ind1, ind2, ind3, o1, o2, knots1, knots2, knots3, pl1, pl2, pl3, uniqtime1, rname1, dtime1, intgtime1, intgname1)
    list(res, knots1, knots2, knots3)
 }


margpartial3 <- function(paras, resp, cov, n, p, cov1, cov2, cov3, Bs1, Bs2, Bs3, A1, A2, A32, A31,ind1, ind2, ind3,   Av1, Av2 ){
    pl <- ncol(A1)
    beta1 <- paras[1 : p]
    beta2 <- paras[(1 + p) : (2 * p)]
    beta3 <- paras[(2 * p + 1) : (3 * p)]
    sp1 <- paras[ (3 * p + 1): (3 * p + pl)]
    sp2 <- paras[ (3 * p + pl + 1): (3 * p + 2 * pl)]
    sp3 <- paras[ (3 * p + 2 * pl + 1): (3 * p + 3 * pl)]
    theta <- paras[  (3 * p + 3 * pl) + 1]
#print(theta)
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
  
    lambda1 <- Bs1 %*% sp1
   
    
    lambda2 <- Bs2 %*% sp2
    
    
    lambda3 <- Bs3 %*% sp3
    
  
    A <- (Lambda1 -  Av1 %*% sp1)  * exp(cov %*% beta1) +  (Lambda2 -  Av2 %*% sp2)  * exp(cov %*% beta2) +  (Lambda32 - Lambda31)  * exp(cov %*% beta3)
    vl <- c( (lambda1 ), (lambda2), (lambda3) )#c(lambda1, lambda2, lambda3)#c(lambda1, lambda2, lambda3)##
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]
#    print(range(theta * A))
#    print(range(vl))
    res <- try((sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* A))+  sumbb + sum(log(vl + 1e-16))))
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-res)
    }
    

    
}






iniestreal3 <- function(theta, bb,  resp, cov){
p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")
      ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
     nml <- max(resp[, 3:5])
     resp[, 3:5] <- resp[, 3:5]
 #resp[, 3:5] <- pnorm(resp[, 3:5])
     knots1 <- quantile(resp[ind1, c(3)], seq(0, 1, length.out = nk))#seq(quantile(resp[ind1, c(3)], 0.2), quantile(resp[ind1, c(3)], 0.8), length.out = nk)
knots1 <- c(unique(max(min(knots1) - 3 * mean(diff(knots1)), 0), max(min(knots1) - 2* mean(diff(knots1)), 0), max(min(knots1) - mean(diff(knots1)), 0)), knots1, max(knots1) + mean(diff(knots1)), max(knots1) + 2 * mean(diff(knots1)), max(knots1) + 3 * mean(diff(knots1)))

     knots2 <-quantile(resp[ind2, c(3)], seq(0, 1, length.out = nk)) #seq(quantile(resp[ind2, c(3)], 0.2), quantile(resp[ind2, c(3)], 0.8), length.out = nk)#quantile(resp[ind2, c(3,  5 )], seq(0.2, 0.8, length = nk))  #
knots2 <- c(unique(max(min(knots2) - 3 * mean(diff(knots2)), 0), max(min(knots2) - 2* mean(diff(knots2)), 0), max(min(knots2) - mean(diff(knots2)), 0)), knots2, max(knots2) + mean(diff(knots2)), max(knots2) + 2 * mean(diff(knots2)), max(knots2) + 3 * mean(diff(knots2)))
     
     knots3 <- quantile(resp[ind3, c(4)], seq(0, 1, length.out = nk))
knots3 <- c(unique(max(min(knots3) - 3 * mean(diff(knots3)), 0), max(min(knots3) - 2* mean(diff(knots3)), 0), max(min(knots3) - mean(diff(knots3)), 0)), knots3, max(knots3) + mean(diff(knots3)), max(knots3) + 2 * mean(diff(knots3)), max(knots3) + 3 * mean(diff(knots3)))#seq(quantile(resp[ind3, c(4)],0.2 ), quantile(resp[ind3, c(4)], 0.8), length.out = nk)#quantile(resp[ind3, c(3,  4)], seq(0.2, 0.8, length = nk)) #quantile(resp[, c(3:5)], seq(0, 1, length = nk))#
uniqtime1 <- unique(as.vector(resp[, c(3, 4, 5)]))
uniqtime1 <- uniqtime1[order(uniqtime1)]
     intervals <- cbind(uniqtime1, c(uniqtime1[-1], 100))[1:(length(uniqtime1)-1), ]
     fill <- function(interval){
         runif(20, interval[1],interval[2])
     }
     intgtime1 <- c(uniqtime1, apply(intervals, 1, fill))
     intgtime1 <- unique(intgtime1)
     intgtime1 <- intgtime1[order(intgtime1)]
     rname1 <- as.character(uniqtime1)
     intgname1 <- as.character(intgtime1)
     dtime1 <- as.matrix(c(0, diff(intgtime1)))
lt <- length(dtime1)
Bs1 <- as.matrix(spline.des(knots1, intgtime1, ord, outer.ok = T)$design)
     Bs2 <-as.matrix( spline.des(knots2, intgtime1, ord, outer.ok = T)$design)
    Bs3 <-as.matrix(spline.des(knots3, intgtime1, ord, outer.ok = T)$design)

   
     
  #   rownames(dA) <- rnA
    mA1 <- apply(diag(as.vector(dtime1)) %*% Bs1, 2, cumsum)
     mA2 <- apply(diag(as.vector(dtime1)) %*% Bs2, 2, cumsum)
mA3 <- apply(diag(as.vector(dtime1)) %*% Bs3, 2, cumsum)
     
     rownames(mA1) <- intgname1
     rownames(mA2) <- intgname1
     rownames(mA3) <- intgname1
      rownames(Bs1) <- intgname1
     rownames(Bs2) <- intgname1
      rownames(Bs3) <- intgname1
     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
    
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
     Av1 <- as.matrix(mA1[as.character(resp[, "v"]), ])
     Av2 <- as.matrix(mA2[as.character(resp[, "v"]), ])
    
     sp1 <- pmax(ginv(t(A1) %*% (A1)) %*% t(A1) %*% bz10$hazard, 0.01)
     sp2 <- pmax(ginv(t(A2) %*% (A2)) %*% t(A2) %*% bz20$hazard, 0.01)
     sp3 <- pmax(ginv(t(A3) %*% (A3)) %*% t(A3) %*% bz30$hazard, 0.01)
     A1 <- as.matrix(mA1[as.character(resp[, "y1"]), ])
     A2 <- as.matrix(mA2[as.character(resp[, "y1"]), ])
     A32 <- as.matrix(mA3[as.character(resp[, "y2"]), ])
     A31 <- as.matrix(mA3[as.character(resp[, "y1"]), ])
Bs1 <- as.matrix(Bs1[as.character(resp[ind1, "y1"]), ])
 Bs2 <- as.matrix(Bs2[as.character(resp[ind2, "y1"]), ])
Bs3 <- as.matrix(Bs3[as.character(resp[ind3, "y2"]), ])


     pl <- ncol(mA1)
     res <- spg(c(beta1, beta2, beta3, sp1, sp2, sp3, theta), margpartial3, gr = NULL, method=2,  lower = c(rep(-1, 3 * p), rep(0, 3 * pl), 0.001), upper = c(rep(100, 3 * p), rep(10, 3 * pl), 10), project=NULL, projectArgs=NULL, control=list(M = 10, maxit = 5000), quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3, Bs1, Bs2, Bs3,  A1, A2, A32, A31, ind1, ind2, ind3,  Av1, Av2)
    res
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
    if(class(res) == "try-error"){
        browser()
    }else{
        return(-res)
    }
    

    
}






iniestreal4 <- function(theta, bb,  resp, cov){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2", "v")
   #  resp[, 3:5] <- exp(resp[, 3:5])/(1 + exp(resp[, 3:5]))
     nml <- max(resp[, 3:5])
     resp[, 3:5] <- resp[, 3:5]
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
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     
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
     ## dA1 <- dA[as.character(resp[, "y1"]), ]
     ## dA2 <- dA[as.character(resp[, "y2"]), ]
     ## dA3 <- dA2
     sp1 <- pmax(ginv(t(A1) %*% (A1)) %*% t(A1) %*% bz10$hazard, 0.01)
     sp2 <- pmax(ginv(t(A2) %*% (A2)) %*% t(A2) %*% bz20$hazard, 0.01)
     sp3 <- pmax(ginv(t(A3) %*% (A3)) %*% t(A3) %*% bz30$hazard, 0.01)
     A1 <- mA1[as.character(resp[, "y1"]), ]
     A2 <- mA2[as.character(resp[, "y1"]), ]
     A32 <- mA3[as.character(resp[, "y2"]), ]
     A31 <- mA3[as.character(resp[, "y1"]), ]
     pl1 <- ncol(mA1)
     pl2 <- ncol(mA2)
     pl3 <- ncol(mA3)
 #   mtheta <<- 0
     
     res <- spg(c(beta1, beta2, beta3, sp1, sp2, sp3, theta), margpartial4, gr = NULL, method=2,  lower = c(rep(-1, 3 * p), rep(0, pl1 + pl2 + pl3), 0.01), upper = c(rep(10, 3 * p), rep(100, pl1 + pl2 + pl3), 10), project=NULL, projectArgs=NULL, control=list(M = 10, maxit = 20000, maxfeval = 1e6, trace = FALSE), quiet=FALSE, resp, cov, n, p, cov1, cov2, cov3, A1, A2, A32, A31, dA1, dA2, dA3, ind1, ind2, ind3,  Av1, Av2)
    res
 }


res <- vector("list")
set.seed(2013)
for(i in 1:100){
    print(i)
    survData <- simCpRsk(1000, p = 1, theta = 0.8, lambda1 = 1, lambda2 = 0.5, lambda3 = 1, kappa = 2,   beta1 = 0.2, beta2 = 0.2, beta3 = 0.2, covm =  NULL, 3, 5)
if(sum(is.na(survData) > 0)){
    next
}
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
res[[i]] <-  iniestreal4(theta, bb, resp, covmy)
 #   print(res[i, 3 * p + 3 * pl + 1])
   #res <- iniestreal2(theta, bb, resp, covmy, knots1, knots2, knots3)
}
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
    t1 <- exp(((u[1]^(-t) - 1) / (t * (lb1 + lb2))))
    t2 <- exp(((u[2]^(-t ) - 1) / (t * (lb1 + lb2))))
    if(r == 0){
        t2 <-  exp((u[3]^(- t / (1 + t)) * ( 1 + t * lb1 * log(t1) + t * lb2 * log(t1)) - ( 1 + t * lb1 * log(t1)  + t * lb2 * log(t1) - t * lb3 * log(t1)))/ (t * lb3))
    }else{
        t1 <- t2 + 3
    }
    c <- runif(1, cen1, cen2) #censoring time
    c(t1, t2, c)
}

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
    c <- runif(1, cen1, cen2)
    #c <- rexp(1, 1/3) #censoring time
    c(t1, t2, c)
}
getres <- function(i){
    res[[i]][[1]]
}
     mres <- do.call(rbind, lapply(1:29, getres))
