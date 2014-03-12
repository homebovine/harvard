
library("survival")
library("BB")
library(numDeriv)
library(rootSolve)
library(parallel)
fA12 <- function(i, beta, vl, resp, cov, n, p){
    ix <- vl[, 2] <= (resp[i, "y1"] + 1e-10)
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(matrix(beta, p, 1)) %*% matrix(cov[i, ], p, 1))
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov, n, p){
    ix1 <- (vl[, 2] > resp[i, "y1"]) & (vl[, 2] <= resp[i, "y2"] + 1e-10)
    
    A3 <- 0
    if(sum(ix1) != 0)
        A3 <- sum(vl[ix1, 1])
    
    A3 <- A3  * exp(t(matrix(beta, p, 1)) %*% matrix(cov[i, ], p, 1))
}
fB <- function(i, theta, resp){
    1/theta + resp[i, "d1"] + resp[i, "d2"]
    
}


getA <- function(i, theta,  beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p){
    
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov, n, p)
    A2 <- fA12(i, beta2, vl2, resp, cov, n, p)
    A3 <- fA3(i, beta3, vl3, resp, cov, n, p)
    A <- A1 + A2 + A3
    Z <- (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + A)
    return(c(A1, A2, A3, A, Z))
}



postz <- function(i, theta, paras, resp){
    (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + paras[i])
}









slvl3 <- function( b1, b2, b3, vl1, vl2, vl3, resp, cov, n, parasall, p, m, f, g, lvl1, lvl2, lvl3){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    zVec <- parasall[5, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec 
    
    fvl <- function(j, lvl, vl,  expb, flg){
        
        a <- (sum(expb[lvl[[j]]]))
        
        if(a != 0 ){
            1/vl[j, 1] - a
        }else{
            browser()
            0
        }
        
        
    }
    
    
    svl1 <- sapply(1 : m, fvl, lvl1, vl1,  expb1, 1) 
    svl2 <- sapply(1 : f, fvl, lvl2, vl2, expb2, 2) 
    svl3 <- sapply(1: g, fvl, lvl3, vl3,  expb3, 0) 
    
    
    
    return(c(svl1, svl2, svl3))
    
}
funcvl <- function(j, vl, flg, resp){
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 1e-10)#&((d1 == 1) )
                                        
        }else if(flg == 2){
             ix <- (y1 >= (vl[j, 2] - 1e-10))# & (d1 == 0)
          
        }else{
           
            ix <- (y2 >= vl[j, 2] - 1e-10) & (y1 < vl[j, 2])#& (d1 == 1)
        }
}
    

warpslv3 <- function(vl, b1, b2, b3, t,  vl1, vl2, vl3, resp, cov, n, svl, parasall, p, m, f, g, lvl1, lvl2, lvl3){
    vl1[, 1] <- (vl[1 : m])
    vl2 [, 1] <- (vl[ (m +1) : (m + f)])
    vl3 [, 1] <- (vl[(m + f + 1) : (m + f + g)])
    paras <- parasall[4, ]
    zVec <- parasall[5, ]
    svl(b1, b2, b3, vl1, vl2, vl3, resp, cov, n, parasall, p, m, f, g, lvl1, lvl2, lvl3)       
}






scoremaxnobb <- function(vbl, resp, cov, beta1, beta2, beta3, theta, vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3){
    vl1<- cbind(exp(vbl[1 : m]), vl1[, 2])
    vl2<- cbind(exp(vbl[(m+1) : (m + f)]), vl2[, 2])
    vl3<- cbind(exp(vbl[(m+f + 1) : (m + f + g)]), vl3[, 2])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl3, parasall, p, m, f, g, lvl1, lvl2, lvl3)
  
    c(vls) 
   
   
}






scoreindv <- function(i, bb, resp, cov, vl1, vl2, vl3, n, p){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
    theta <- exp(bb[3*p + 1])
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov, n, p)
    A2 <- fA12(i, beta2, vl2, resp, cov, n, p)
    A3 <- fA3(i, beta3, vl3, resp, cov, n, p)
    A <- A1 + A2 + A3
    x <- cov[i, ]
    B <- fB(i, theta, resp)
    commd <- (1 + theta * A)
    if(is.na(log(commd))){
        browser()
        print(c(A, theta))
    }
    u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
    u2 <- d1 * x - B * x * theta * A1 / commd ##a p \times 1 vector
    u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
    u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
    mu <- c(u2, u3, u4, u1) #(3 * p + 1)vector
    return(mu)
    


}

scorepara <- function(i, bb, resp, cov, vl1, vl2, vl3, n, p){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
    theta <- exp(bb[3*p + 1])
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    x <- cov[i, ]
    A1 <- exp( sum(beta1 * x) ) * resp[i, "y1"] ^ 2
    A2 <- exp( sum(beta2 * x)) * resp[i, "y1"] ^ 2
    A3 <- exp( sum(beta3 * x)) * (resp[i, "y2"] ^ 2 - resp[i, "y1"]^2)
    A <- A1 + A2 + A3
    
    B <- fB(i, theta, resp)
    commd <- (1 + theta * A)
    if(is.na(log(commd))){
        browser()
        print(c(A, theta))
    }
    u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
    u2 <- d1 * x - B * x * theta * A1 / commd ##a p \times 1 vector
    u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
    u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
    mu <- c(u2, u3, u4, u1) #(3 * p + 1)vector
    return(mu)
    


}

wrapscoreindv <- function(bb, resp, cov, n, vl1, vl2, vl3,  p){
    apply(sapply(1 : n, scoreindv, bb, resp, cov, vl2, vl2, vl3, n, p), 1, sum)[1 : (3 * p+1)]
    }
wrapscoreindv1 <- function(bb, resp, cov, n, vl1, vl2, vl3,  p){
    apply(sapply(1 : n, scorepara, bb, resp, cov, vl2, vl2, vl3, n, p), 1, sum)
    }



    
margpartial <- function(vlbb, vl1, vl2, vl3, resp, cov, n, p ){
    vl <- vlbb[1 : (m + f + g)]
    vl1[, 1] <- vlbb[1 :m ]
    vl2[, 1] <- vlbb[(m + 1) : (m + f) ]
    vl3[, 1] <- vlbb[(m + f + 1) : (m + f + g) ]
    beta1 <- vlbb[(m + f + g + 1) : (m + f + g + p) ]
    beta2 <- vlbb[(m + f + g + p + 1) : (m + f + g + 2* p) ]
    beta3 <- vlbb[(m + f + g + 2 * p + 1) : (m + f + g + 3 * p)]
    theta <- exp(vlbb[m + f + g + 3 * p + 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    paras <- parasall[, 4]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    B <- 1/theta + resp[, 1] + resp[, 2]  
    sum(B1) * log(theta + 1)  - sum(B * log(1 + theta* paras))+  sumbb + sum(log(vl))
}



getfromlist<- function(itr, lst, ind){
    lst[[itr]][[ind]]
}
getvlist <- function(i, estres, ix, ix1){
    estres[[ix]][[2]][[ix1]][, 1]
}
comscore <- function(vlbb, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3){
    vl <- (vlbb[1 : (m + f + g)])
    vl1[, 1] <- exp(vlbb[1 :m ])
    vl2[, 1] <- exp(vlbb[(m + 1) : (m + f) ])
    vl3[, 1] <- exp(vlbb[(m + f + 1) : (m + f + g) ])
    beta1 <- vlbb[(m + f + g + 1) : (m + f + g + p) ]
    beta2 <- vlbb[(m + f + g + p + 1) : (m + f + g + 2* p) ]
    beta3 <- vlbb[(m + f + g + 2 * p + 1) : (m + f + g + 3 * p)]
    ltheta <- (vlbb[m + f + g + 3 * p + 1])
    nobb <- scoremaxnobb(vl, resp, cov, beta1, beta2, beta3, exp(ltheta), vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3) 
    novl <- wrapscoreindv(c(beta1, beta2, beta3, (ltheta)), resp, cov, n, vl1, vl2, vl3, p) 
    c(nobb, novl) /n
    
}
comscore1 <- function(vlbb, vl1, vl2, vl3, resp, cov,  n, p, m, f, g){
    vl <- (vlbb[1 : (m + f + g)])
    vl1[, 1] <- exp(vlbb[1 :m ])
    vl2[, 1] <- exp(vlbb[(m + 1) : (m + f) ])
    vl3[, 1] <- exp(vlbb[(m + f + 1) : (m + f + g) ])
    beta1 <- vlbb[(m + f + g + 1) : (m + f + g + p) ]
    beta2 <- vlbb[(m + f + g + p + 1) : (m + f + g + 2* p) ]
    beta3 <- vlbb[(m + f + g + 2 * p + 1) : (m + f + g + 3 * p)]
    ltheta <- (vlbb[m + f + g + 3 * p + 1])
    nobb <- scoremaxnobb(vl, resp, cov, beta1, beta2, beta3, exp(ltheta), vl1, vl2, vl3, n, p, m, f, g) 
    novl <- wrapscoreindv(c(beta1, beta2, beta3, (ltheta)), resp, cov, n, vl1, vl2, vl3, p)
    c(nobb, novl) 
    
}
simCpRsk <- function(n, p, l1, l2, l3, a, b1, b2, b3, c1, c2, cen1, cen2){
    covm <- matrix(runif(p * n, c1, c2), n, p) #covariance matrix 
    simdata1 <- t(sapply(1 : n, simwei2, theta,  l1, l2, l3, b1, b2, b3, a, covm, cen1, cen2))
    d1 <- simdata1[, 1] < simdata1[, 2]&simdata1[, 1] < simdata1[, 3]
    d2 <- simdata1[, 2] < simdata1[, 3]
    y2 <- pmin(simdata1[, 2], simdata1[, 3])
    y1 <- pmin(simdata1[, 1], y2)
    simresp1 <- cbind(d1, d2, y1, y2)
    nsimresp1 <- simdata1
    colnames(simresp1) <- c("d1", "d2", "y1", "y2")
    return(simresp1, covm)
}
simeval <- function(sitr){
   res <-  estreal(c(0, b1, 0, b2, 0, b3, log(2)), lsimresp1[[sitr]], lcovm[[sitr]], FALSE, 1, 1e-3, 1e-3)[[1]]
}
