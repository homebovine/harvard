
                                        #vl1 is a m x 2 dimensional vector, vl1[, "lambda value"], vl1[, "time"]
library("survival")
library("BB")
library(numDeriv)
library(rootSolve)
fA12 <- function(i, beta, vl, resp, cov){
    ix <- vl[, 2] <= (resp[i, "y1"] + 1e-10)
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(matrix(beta, p, 1)) %*% matrix(cov[i, ], p, 1))
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov){
    ix1 <- (vl[, 2] > resp[i, "y1"]) & (vl[, 2] <= resp[i, "y2"] + 1e-10)
    
    A3 <- 0
    if(sum(ix1) != 0)
        A3 <- sum(vl[ix1, 1])
    
    A3 <- A3  * exp(t(matrix(beta, p, 1)) %*% matrix(cov[i, ], p, 1))
}
fB <- function(i, theta, resp){
    1/theta + resp[i, "d1"] + resp[i, "d2"]
    
}


getA <- function(i, theta,  beta1, beta2, beta3, vl1, vl2, vl3){
    
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov)
    A2 <- fA12(i, beta2, vl2, resp, cov)
    A3 <- fA3(i, beta3, vl3, resp, cov)
    A <- A1 + A2 + A3
    Z <- (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + A)
    return(c(A1, A2, A3, A, Z))
}



postz <- function(i, theta, paras, resp){
    (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + paras[i])
}




slvl31 <- function( b1, b2, b3, vl1, vl2, vl3, resp, cov, n, parasall){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    zVec <- parasall[5, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec 
    
    fvl <- function(j, vl,  expb, flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 1e-10)
            if(0){
            b <- sum(y1 <= (vl[j, 2] + 1e-10) & y1 >  (vl[j - 1, 2] ))
            }else{
                b <- 1
                }
                                        
        }else if(flg == 2){
             ix <- (y1 >= (vl[j, 2] - 1e-10))# & (d1 == 0)
             if(0){
             b <- sum((y2 <= (vl[j, 2]) + 1e-10) & (d1 == 0)&y2 >  (vl[j - 1, 2] ))
             }else{
                 b <- 1
                 }
        }else{
           
            ix <- (y2 >= vl[j, 2] - 1e-10) & (y1 < vl[j, 2])#& (d1 == 1)
            if(0){
                 b <- sum((y2 <= (vl[j, 2]) + 1e-10) & (d1 == 1)&y2 >  (vl[j - 1, 2] ))
                }else{
                    b <- 1
                   }
        }
        
        a <- (sum(expb[ix]))
        
        if(a != 0 ){
            b/ a
        }else{
            browser()
            0
        }
        
        
    }    
    
    svl1 <- sapply(1 : m, fvl, vl1, expb1, 1)
    svl2 <- sapply(1 : f, fvl, vl2, expb2, 2)
    svl3 <- sapply(1: g, fvl, vl3,  expb3, 0)
    
    
    
    return(c(svl1, svl2, svl3))
    
}



slvlb <- function( bb,   vl1, vl2, vl3, resp, cov, n, parasall){
    b1 <- matrix(bb[1 :p], p, 1)
    b2 <- matrix(bb[(p + 1):(2*p)], p, 1)
    b3 <- matrix(bb[(2 * p+ 1):(3*p)], p, 1)
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    paras <- parasall[4, ]
    zVec <- parasall[5, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec
    expzb1 <- matrix(expb1, n, p) * cov
    expzb2 <- matrix(expb2, n, p) * cov
    expzb3 <- matrix(expb3, n, p) * cov
    
    fvl <- function(j, vl,  expb, expzb, flg, cov){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 1e-10)
        }else 
            ix <- (y2 >= vl[j, 2] - 1e-10) & (y1 < vl[j, 2])
        
        a <- (sum(expb[ix]))
        mxe <- matrix(expzb[ix, ], ncol = p)
        b <- apply(mxe, 2, sum)
        
        if(a != 0 ){
            cov[j, ] - b /a
        }else{
            browser()
            0
        }
        
        
    }
    
    
    svl1 <- apply(matrix(sapply(1 : m, fvl, vl1, expb1, expzb1, 1, cov1), nrow = p), 1, sum)
    svl2 <- apply(matrix(sapply(1 : f, fvl, vl2, expb2, expzb2, 1, cov2), nrow = p), 1, sum)
    svl3 <- apply(matrix(sapply(1: g, fvl, vl3,  expb3, expzb3, 0, cov3), nrow = p), 1, sum)
    
    
    
    return(c(svl1, svl2, svl3))
    
}

slvl3 <- function( b1, b2, b3, vl1, vl2, vl3, resp, cov, n, parasall){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    zVec <- parasall[5, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec 
    
    fvl <- function(j, vl,  expb, flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 1e-10)#&((d1 == 1) )
          #  if(0){
           # b <- sum(y1 <= (vl[j, 2] + 1e-10) & y1 >  (vl[j - 1, 2]))
           # }else{
            #    b <- 1
             #   }
                                        
        }else if(flg == 2){
             ix <- (y1 >= (vl[j, 2] - 1e-10))# & (d1 == 0)
          #   if(0){
           #  b <- sum((y2 <= (vl[j, 2]) + 1e-10) & (d1 == 0)&y2 >  (vl[j - 1, 2] ))
            # }else{
             #    b <- 1
              #   }
        }else{
           
            ix <- (y2 >= vl[j, 2] - 1e-10) & (y1 < vl[j, 2])#& (d1 == 1)
            ## if(0){
            ##      b <- sum((y2 <= (vl[j, 2]) + 1e-10) & (d1 == 1)&y2 >  (vl[j - 1, 2]))
            ##     }else{
            ##         b <- 1
       #            }
        }
        
        a <- (sum(expb[ix]))
        
        if(a != 0 ){
            1/vl[j, 1] - a
        }else{
            browser()
            0
        }
        
        
    }
    
    
    svl1 <- sapply(1 : m, fvl, vl1, expb1, 1)
    svl2 <- sapply(1 : f, fvl, vl2, expb2, 2)
    svl3 <- sapply(1: g, fvl, vl3,  expb3, 0)
    
    
    
    return(c(svl1, svl2, svl3))
    
}


warpslv3 <- function(vl, b1, b2, b3, t,  vl1, vl2, vl3, resp, cov, n, svl, parasall){
    vl1[, 1] <- (vl[1 : m])
    vl2 [, 1] <- (vl[ (m +1) : (m + f)])
    vl3 [, 1] <- (vl[(m + f + 1) : (m + f + g)])
    paras <- parasall[4, ]
    zVec <- parasall[5, ]
    svl(b1, b2, b3, vl1, vl2, vl3, resp, cov, n, parasall)   
    
    }


slvlb2 <- function( bb,   vl1, vl2, vl3, resp, cov, n, parasall){
    b1 <- matrix(bb[1 :p], p, 1)
    b2 <- matrix(bb[(p + 1):(2*p)], p, 1)
    b3 <- matrix(bb[(2 * p+ 1):(3*p)], p, 1)
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    paras1 <- parasall[1, ]
    paras2 <- parasall[2, ]
    paras3 <- parasall[3, ]
    zVec <- parasall[5, ]
    expb1 <- paras1  * zVec
    expb2 <- paras2  * zVec
    expb3 <- paras3  * zVec
    expzb1 <- matrix(expb1, n, p) * cov
    expzb2 <- matrix(expb2, n, p) * cov
    expzb3 <- matrix(expb3, n, p) * cov
    
    svl1 <- t(d1) %*% (cov ) - apply(expzb1, 2, sum)
    svl2 <- t((1-d1) * d2) %*% (cov) - apply(expzb2, 2, sum)
    svl3 <- t((d1) * d2) %*% (cov) - apply(expzb3, 2, sum)
    return(c(svl1, svl2, svl3))
    
}
wrapslvlb2 <- function(bb, vl1, vl2, vl3, resp, cov, n, parasall){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2 * p + 1):(3 * p)]
    t <- (bb[3 * p + 1])#exp(theta)
    #bb <- c(beta1, beta2, beta3, t)
    paras <- parasall[4, ]
    zVec <- parasall[5, ]
    if(is.na(log(t))){
        browser()
        }
    logZ <- sapply(1:n, postlogz, t, paras, resp)
   dilogZ <- sum(logZ) - sum(zVec)
    
    
    sb <- slvlb(bb[1 : (3*p)], vl1, vl2, vl3, resp, cov, n, parasall)
    st <- slvtheta(t, dilogZ)
    as.numeric(c(sb, st))
}



scoremax <- function(vbl){
    vl1[, 1]<- (vbl[1 : m])
    vl2[, 1]<- (vbl[(m+1) : (m + f)])
    vl3[, 1]<- (vbl[(m+f + 1) : (m + f + g)])
  #  nvl1 <- vl1
  #  nvl2 <- vl2
  #  nvl3 <- vl3
    beta1 <- vbl[(m + f+ g+ 1): (m +f+ g+ p)]
    beta2 <- vbl[(m + f+ g+ p + 1): (m + f+ g+ 2* p)]
    beta3 <- vbl[(m + f+ g+ 2* p + 1): (m + f+ g+ 3* p)]
   # nvl1[, 1]<- vl1[, 1] * exp(beta1 * cov1)
   # nvl2[, 1]<- vl2[, 1] * exp(beta2 * cov2)
   # nvl3[, 1]<- vl3[, 1] * exp(beta3 * cov3)
    theta <- (vbl[length(vbl)])
   
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl3, parasall)
    bsc <-  wrapslvlb2(c(beta1, beta2, beta3, (theta)), vl1, vl2, vl3, resp, cov, n, parasall)
    -c(vls, bsc)
}
scoremax1 <- function(vbl){
    #vl1[, 1]<- (vbl[1 : m])
    #vl2[, 1]<- (vbl[(m+1) : (m + f)])
    #vl3[, 1]<- (vbl[(m+f + 1) : (m + f + g)])
    beta1 <- vbl[(1): (p)]
    beta2 <- vbl[( p+ 1): (2* p)]
    beta3 <- vbl[(2* p + 1): (3* p)]
    theta <- (vbl[3*p + 1]) + 1e-10
   # print(theta)
    crit <- 0
    zVec <- 0
    for(i in 1: 200){
        #print(c(crit, i))
        parasall <- sapply(1 : n, getA, theta,   0, 0, 0, nvl1, nvl2, nvl3)
        crit <- sum((zVec - parasall[5, ])^2)
        vls <- warpslv3(c((nvl1[, 1]), (nvl2[, 1]), (nvl3[, 1])), 0, 0, 0, theta, nvl1, nvl2, nvl3, resp, cov, n, slvl31, parasall)
        nvl1[, 1]<- (vls[1 : m])
        nvl2[, 1]<- (vls[(m+1) : (m + f)])
        nvl3[, 1]<- (vls[(m+f + 1) : (m + f + g)])
        zVec <- parasall[5, ]
        if(crit <= 1e-7){
            break
            }
                
                
        
        }
    
    #parasall <- sapply(1 : n, getA, theta,   0, 0, 0, vl1, vl2, vl3)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl31, parasall)
    vl1[, 1]<- (vls[1 : m])
    vl2[, 1]<- (vls[(m+1) : (m + f)])
    vl3[, 1]<- (vls[(m+f + 1) : (m + f + g)])
    print(theta)
    bsc <-  wrapslvlb2(c(beta1, beta2, beta3, theta), vl1, vl2, vl3, resp, cov, n, parasall)
   # -c(vls, bsc)
}
wraplike <- function(vbl){
    vl1[, 1]<- (vbl[1 : m])
    vl2[, 1]<- (vbl[(m+1) : (m + f)])
    vl3[, 1]<- (vbl[(m+f + 1) : (m + f + g)])
    nvl1 <- vl1
    nvl2 <- vl2
    nvl3 <- vl3
     beta1 <- vbl[(m + f+ g+ 1): (m +f+ g+ p)]
    beta2 <- vbl[(m + f+ g+ p + 1): (m + f+ g+ 2* p)]
    beta3 <- vbl[(m + f+ g+ 2* p + 1): (m + f+ g+ 3* p)]
    
    nvl1[, 1]<- vl1[, 1] * exp(beta1 * cov1)
    nvl2[, 1]<- vl2[, 1] * exp(beta2 * cov2)
    nvl3[, 1]<- vl3[, 1] * exp(beta3 * cov3)
    theta <- (vbl[length(vbl)])
   
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
    paras <- parasall[4, ]
   # vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl31, parasall)
   # vl1[, 1]<- vls[1 : m]
   # vl2[, 1]<- vls[(m+1) : (m + f)]
    #vl3[, 1]<- vls[(m+f + 1) : (m + f + g)]
    
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}
wraplike1 <- function(vbl){
    #vl1[, 1]<- vbl[1 : m]
    #vl2[, 1]<- vbl[(m+1) : (m + f)]
    #vl3[, 1]<- vbl[(m+f + 1) : (m + f + g)]
    beta1 <- vbl[( 1): ( p)]
    beta2 <- vbl[( 1): ( 2* p)]
    beta3 <- vbl[(2* p + 1): (3* p)]
    theta <- vbl[length(vbl)]
    parasall <- sapply(1 : n, getA, theta,  0, 0, 0, vl1, vl2, vl3)
    paras <- parasall[4, ]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}
wraplikenovl <- function(vbl){
    beta1 <- vbl[( 1): ( p)]
    beta2 <- vbl[( p + 1): ( 2* p)]
    beta3 <- vbl[( 2* p + 1): (3* p)]
 #   theta <- vbl[length(vbl)]
    paras <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)[4, ]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}

scoremaxnovl <- function(vbl){
    beta1 <- vbl[( 1): (p)]
    beta2 <- vbl[( p + 1): ( 2* p)]
    beta3 <- vbl[(2* p + 1): (3* p)]
    theta <- exp(vbl[3 * p + 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
     bsc <-  wrapslvlb2(c(beta1, beta2, beta3, theta), vl1, vl2, vl3, resp, cov, n, parasall)
  #  bsc <-  wrapslvlb2(c(beta1, beta2, beta3, log(theta)), vl1, vl2, vl3, resp, cov, n)
    -c( bsc)
   
   
}

scoremaxnobb <- function(vbl){
    vl1[, 1]<- exp(vbl[1 : m])
    vl2[, 1]<- exp(vbl[(m+1) : (m + f)])
    vl3[, 1]<- exp(vbl[(m+f + 1) : (m + f + g)])
    #vl1[, 1]<- approxfun(nvl1[, 2], nvl1[, 1])(vl1[, 2])
    #vl2[, 1]<- approxfun(nvl2[, 2], nvl2[, 1])(vl2[, 2])
    #vl3[, 1]<- approxfun(nvl3[, 2], nvl3[, 1])(vl3[, 2])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl3, parasall)
   # bsc <-  wrapslvlb2(c(beta1, beta2, beta3, log(theta)), vl1, vl2, vl3, resp, cov, n)
    -c(vls)
   
   
}
wraplikenobb <- function(vbl){
    vl1[, 1]<- vbl[1 : m]
    vl2[, 1]<- vbl[(m+1) : (m + f)]
    vl3[, 1]<- vbl[(m+f + 1) : (m + f + g)]
    paras <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)[4, ]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}

scoreindv <- function(i, bb, resp, cov){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
    theta <- exp(bb[3*p + 1])
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov)
    A2 <- fA12(i, beta2, vl2, resp, cov)
    A3 <- fA3(i, beta3, vl3, resp, cov)
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

wrapscoreindv <- function(bb, resp, cov, n){
    apply(sapply(1 : n, scoreindv, bb, resp, cov), 1, mean)
    }


    
margpartial <- function(theta, vl, paras, sumbb ){
    B <- 1/theta + B2    
    sum(B1) * log(theta + 1)  - sum(B * log(1 + theta* paras))+  sumbb + sum(log(vl))
}


slvtheU <- function(theta, paras, resp){
    sum( resp[, 1] * resp[, 2] / ( (theta + 1)) + 1/(theta^2) * log(1 + theta * paras) - (1 / theta + resp[, 1] + resp[, 2]) * paras / ( 1 + theta * paras))
}

slvtheta <- function(t, dilogZ){
    dilogZ - n* ( digamma(1/t) ) - n * (log(t) - 1)
}
postlogz <- function(i, theta, paras, resp){
    digamma(1/theta + resp[i, 1]+ resp[i, 2]) - log(1/theta + paras[i])
}





mtbb <- 1
nitr <- 10
mres <- mvl <- vector("list")

for(simitr in 1 : nsim){
    theta <- 1
    print(simitr)
    cov <- matrix(lcovm[[simitr]], n, p)
    resp <- lsimresp1[[simitr]]
    B1 <- resp[, 1] * resp[, 2]
    B2 <- resp[, 1] + resp[, 2]   
    colnames(resp) <- c("d1", "d2", "y1", "y2")
    d1 <- resp[, 1]
    d2 <- resp[, 2]
    y1 <- resp[, 3]
    y2 <- resp[, 4]
    ind1 <- which(d1 == 1 )
    ind2 <- which((d1 == 0) & (d2 == 1))
    ind3 <- which((d1 ==1) & (d2 == 1))
    p <- ncol(cov)
    subgix1 <- (d1 == 0)
    respsub1 <- resp[subgix1, ]
    respsub2 <- resp[ind1, ]

    surv1 <- (coxph(Surv(resp[, "y1"], resp[, "d1"]) ~ (cov)))
    surv2 <- (coxph(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~ cov[subgix1]))
    surv3 <- (coxph(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~ cov[ind1, ]))
    o1 <- order(resp[ind1, 3])
    o2 <- order(resp[ind2, 4])
    o3 <- order(resp[ind3, 4])
    cov1 <- cov[ind1, , drop = F][o1, , drop = F]
    cov2 <- cov[ind2, , drop = F][o2, , drop = F]
    cov3 <- cov[ind3, , drop = F][o3, , drop = F]
    m0 <- length(ind1)
    f0 <- length(ind2)
    g0 <- length(ind3)
    m <- m0
    f <- f0
    g <- g0
    ## q1 <- quantile(resp[ind1, 3], seq(0, 1, length.out = m), type = 3)
    ## q2 <- quantile(resp[ind2, 4], seq(0, 1, length.out = f), type = 3)
    ## q3 <- quantile(resp[ind3, 4], seq(0, 1, length.out = g), type = 3)
    ## sbm <- which(resp[ind1, 3][o1] %in% q1)#sample(1 : m, 10)
    ## sbf <- which(resp[ind2, 4][o2] %in% q2)
    ## sbg <- which(resp[ind3, 4][o3] %in% q3)
    ## sbm <- sbm[order(sbm)]
    ## sbf <- sbf[order(sbf)]
    ## sbg <- sbg[order(sbg)]
    ## cov1<- cov1[sbm, , drop = F]
    ## cov2 <- cov2[sbf, , drop = F]
    ## cov3 <- cov3[sbg, , drop = F]
    
    vl10 <- matrix(NA, m, 2)
    vl20 <- matrix(NA, f, 2)
    vl30 <- matrix(NA, g, 2)
    bz10 <- basehaz(surv1)
    bz20 <- basehaz(surv2)
    bz30 <- basehaz(surv3)
    ixbz1 <- which(c(bz10[1, 1], diff(bz10[, 1])) != 0)
    ixbz2 <- which(c(bz20[1, 1], diff(bz20[, 1])) != 0)
    ixbz3 <- which(c(bz30[1, 1], diff(bz30[, 1])) != 0)
    bz1 <- bz10[ixbz1, ]#[sbm, ]
    bz2 <- bz20[ixbz2, ]#[sbf, ]
    bz3 <- bz30[ixbz3, ]#[sbg, ]
    bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
    bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
    bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))
    
    vl10 <- bz1
    vl20 <- bz2
    vl30 <- bz3
    bz1 <- bz10[ixbz1, ][sbm, ]
    bz2 <- bz20[ixbz2, ][sbf, ]
    bz3 <- bz30[ixbz3, ][sbg, ]
    bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
    bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
    bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))
    nvl10 <- bz1
    nvl20 <- bz2
    nvl30 <- bz3
    m <- m0#nrow(nvl10)
    f <- f0#nrow(nvl20)
    g <- g0#nrow(nvl30)
    
    ## theta <- 1
    ## beta1 <- surv1$coefficients
    ## beta2 <- surv2$coefficients
    ## beta3 <- surv3$coefficients
#    bb <- c(beta1, beta2, beta3, theta)
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
    vl3 <- vl30
    nvl1 <- nvl10
    nvl2 <- nvl20
    nvl3 <- nvl30
    vl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]))
    beta1 <- 1
    beta2 <- 1
    beta3 <- 0.5
    theta <- 1
    vl1[, 1] <- vl[1 : m]
    vl2[, 1] <- vl[ (m + 1): (m + f)]
    vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
    # vl1[, 1]<- approxfun(nvl1[, 2], nvl1[, 1])(vl1[, 2])
    # vl2[, 1]<- approxfun(nvl2[, 2], nvl2[, 1])(vl2[, 2])
    # vl3[, 1]<- approxfun(nvl3[, 2], nvl3[, 1])(vl3[, 2])
     parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
     paras <- parasall[4, ]
     crit <- 0
     broot0 <- broot <- c(0, 0, 0, -0.5)
     for(i in 1){
         evl <- dfsane(c(log(vl10[, 1]), log(vl20[, 1]), log(vl30[, 1])), scoremaxnobb, method = 2, control = list(trace = FALSE), quiet = FALSE)$par
         vl <- exp(evl)
         vl1[, 1] <- vl[1 : m]
         vl2[, 1] <- vl[ (m + 1): (m + f)]
         vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
         brootn<- multiroot(wrapscoreindv, broot, maxiter = 100, rtol = 1e-06, atol = 1e-08, ctol = 1e-08, useFortran = TRUE, positive = FALSE, jacfunc = NULL,  jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,  resp, cov, n)$root
         crit <- sum((brootn /4- broot/4)^2)
         
         if(crit <= 1e-4){
             break
         }else{
             broot <- brootn
             beta1 <- brootn[1]
             beta2 <- brootn[2]
             beta3 <- brootn[3]
             theta <- exp(brootn[4])
              print(crit)
         }
             
                                       
        
}
   # temp <- spg(c(vl1[, 1], vl2[, 1], vl3[, 1], beta1, beta2, beta3, (theta)), wraplike, gr = scoremax, method = 3, project = NULL, lower = c(rep(0.0001, m + f+ g), rep(-10, 3 * p), 0.0001), upper = c(rep(2, m + f+ g), rep(5, 3 * p), 2), projectArgs = NULL, control = list( ftol = 1e-5, gtol = 1.e-03, maxit = 1500, checkGrad = FALSE), quiet = FALSE)
    

      #                                 temp <-  multiroot(scoremax1, c( beta1, beta2, beta3, 1), maxiter = 100, rtol = 1e-06, atol = 1e-05, ctol = 1e-05, useFortran = TRUE, positive = TRUE, jacfunc = NULL, jactype = "fullint",  verbose = FALSE, bandup = 1, banddown = 1)  #sane(c(log(vl1[, 1]), log(vl2[, 1]), log(vl3[, 1]),  beta1, beta2, beta3, log(theta)), scoremax, method = 2, control = list(), quiet = FALSE) #
    
    #vbl <- temp$par
    #vl1[, 1]<- (vbl[1 : m])
   # vl2[, 1]<- (vbl[(m+1) : (m + f)])
    #vl3[, 1]<- (vbl[(m+ f + 1) : (m + f + g)])
    ## beta1 <- vbl[(m + f + g + 1): (m + f + g+ p)]
    ## beta2 <- vbl[(m + f+ g + p + 1): (m + f + g + 2* p)]
    ## beta3 <- vbl[(m + f+ g + 2* p + 1): (m + f + g + 3* p)]
    ## theta <- (vbl[length(vbl)])

    print(c(beta1, beta2, beta3, theta))

    mres[[simitr]] <- c(beta1, beta2, beta3, theta, crit)
    mvl[[simitr]] <- list(vl1, vl2, vl3)

}

############simulate samples
integrand <- function(t2){
    a2 * l2 * t2 ^ (a2 - 1) * (1 + t * l1 * t2 ^ a1 + t * l2 * t2 ^ a2)^(-1/t - 1)
    }
uncd2 <- function(t2){
    integrate(integrand, t2, Inf)$value
    }
  p <-   integrate( integrand, 0, Inf)$value
integrand1 <- function(t2, t1){
    (t + 1) * a1 * l1 * t1 ^ (a1 - 1) * a2 * l2 * t2 ^ (a2 - 1) * (1 + t * l1 * t2 ^ a1 + t * l3* t2 ^ a3)^(-1/t - 2)
    }
cd21 <- function(t2, t1){
   integrate(integrand1, max(t2, t1), Inf)$value/integrate(integrand1, t1, Inf)$value
 }
cd1 <- function(t1){
    integrate(integrand1, t1, Inf)$value
    
    }


simufun <- function(i, t, a1, a2, a3, l1, l2, l3, b1, b2, b3, cov){
g <- rgamma(1, shape = 1/t, rate = 1/t)
r <- rbinom(1, 1, 0.5)
u <- -log(runif(3))
#u3 <- runif(1)
expb1 <- exp(sum(cov[i, ] * b1))
expb2 <- exp(sum(cov[i, ] * b2))
expb3 <- exp(sum(cov[i, ] * b3))
lb1 <- g *l1 * expb1
lb2 <- g * l2 * expb2
lb3 <- g * l3 * expb3
t1 <- (u[1]/ (lb1))^(1/a1)
t2 <- (u[2]/ (lb2))^(1/a2)
integrand <- function(t1, t){
    exp(- lb3 * (t + t1)^ a3) * a1 * lb1 * t1 ^ (a1 - 1) * exp(-lb1 * t1 ^ a1)
    }
rootfun <- function(t){
    integrate( integrand, 0, Inf, t)$value - u3
    }

if(t1 < t2){
    #t2 <- t1 + uniroot(rootfun, c(0, 10000))$root
  #  t2 <-   ((u[3] + lb3 * t1 ^ a3)/ ((lb3)))^(1/a3)
    t2 <-   t1 + ((u[3])/ (lb3))^(1/a3)
}else{
    t1 <- t2 + 3
}

c <- runif(1, 2.5, 5)
c(t1, t2, c)
}

simwei <- function(i, t, l1, l2, l3, b1, b2, b3, a, cov){
    
    
    expb1 <- exp(sum(cov[i, ] * b1))
    expb2 <- exp(sum(cov[i, ] * b2))
    expb3 <- exp(sum(cov[i, ] * b3))
    lb1 <- l1 * expb1
    lb2 <- l2 * expb2
    lb3 <- l3 * expb3
    p <- lb2/(lb1 + lb2)
    r <- rbinom(1, 1, p)
    u <- runif(3) 
    t1 <- ((u[1]^(-(t *(lb1 + lb2)) /lb1 ) - 1) / (t * (lb1 + lb2)))^(1/a)
    t2 <- ((u[2]^(-(t *(lb1 + lb2)) /lb2 ) - 1) / (t * (lb1 + lb2)))^(1/a)
    if(t1 < t2 ){
        t1 <- ((u[1]^(-(t *(lb1 + lb2)) /lb1 ) - 1) / (t * (lb1 + lb2)))^(1/a)
        t2 <-  ((u[3]^(- t / (1 + t)) * ( 1 + t * lb1 * t1 ^ a + t * lb2 * t1^a) - ( 1 + t * lb1 * t1 ^ a + t * lb2 * t1^a - t * lb3 * t1 ^ a))/ (t * lb3))^ (1/a)
    }else{
        
        t1 <- t2 + 3
        }
    c <- runif(1, 2.5, 5)
    c(t1, t2, c)
    }
simwei2 <- function(i, t, l1, l2, l3, b1, b2, b3, a, cov){
    
    
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
    c <- runif(1, 2.5, 5)
    c(t1, t2, c)
    }

simunonhom <- function(i, t, l1, l2, l3, b1, b2, b3, cov){
    g <- rgamma(1, shape = 1/t, rate = 1/t)
    mb <- (b1 + b2 + b3) / 3
    pp <- 1# exp(sum(cov[i, ] * mb)) 
    r <- rbinom(1, 1, pp/(1 + pp))
    expb1 <- 2 * g * exp(sum(cov[i, ] * b1))
    expb2 <- 2 * g * exp(sum(cov[i, ] * b2))
    expb3 <- 4 * g * exp(sum(cov[i, ] * b3))
    u <- runif(3)
    ch <- -log(u)
    ch[1] <- ch[1]/ (expb1)
    ch[2] <- ch[2]/ (expb2)
    if(ch[1]> 1- exp(-3) ){
        t1 <- (ch[1]- (1 - exp(-3))) /exp(-3) + 3
    }else
        t1 <- (-log(1 - ch[1]))
    if(ch[2]> 1- exp(-3) ){
        t2 <- (ch[2] - (1 - exp(-3))) /exp(-3) + 3
    }else
        t2 <- (-log(1 - ch[2]))
    

    if(r == 1){
        if(t1 < 3){
            s <- 0#(expb3 * (1- exp(-t1)))
            chs <- (ch[3] +s)/( expb3)
            if(chs  > exp(-t1)- exp(-3) ){
                t2 <- (chs - (exp(-t1) - exp(-3))) /exp(-3) + 3
            }else
                t2 <- (-log(exp(-t1) - chs))
            
        }else{
            s <- 0#(expb3 * (1- exp(-3) + exp(-3)*(t1 -3)))
            chs <- (ch[3] +s)/(expb3)
            t2 <- chs /(exp(-3)) + t1# (chs - (exp(-t1) - exp(-3))) /exp(-3) + 3#
        }
    
    
    }else{
        t1 <- t2 + 3
    }
    r <- rbinom(1, 1, 0.5)
    if(r == 1)
        c <- runif(1, 1.5, 3)
    else
        c<- 3
    c(t1, t2, c)
}

l1 <- 1
l2 <- 1
l3 <- 2
a1 <- 2
a2 <- 2
a3 <- 2
b1 <- c(1)
b2 <- c(1)
b3 <- c(0.5)

lsimresp1 <- lsimnresp1 <- lcovm <-  vector("list")
set.seed(2014)
for(itr in 1 : nsim){
    covm <- matrix(runif(p * n, 0, 0.5), n, p)
    simdata1 <- t(sapply(1 : n, simwei2, 1,  l1, l2, l3, b1, b2, b3, a1, covm))
    d1 <- simdata1[, 1] < simdata1[, 2]&simdata1[, 1] < simdata1[, 3]
    d2 <- simdata1[, 2] < simdata1[, 3]
    y2 <- pmin(simdata1[, 2], simdata1[, 3])
    y1 <- pmin(simdata1[, 1], y2)
    simresp1 <- cbind(d1, d2, y1, y2)
    nsimresp1 <- simdata1
    colnames(simresp1) <- c("d1", "d2", "y1", "y2")

    lsimresp1[[itr]] <- simresp1
    lsimnresp1[[itr]] <- simdata1
    lcovm[[itr]] <- covm
}




vcr <- 0
for(theta in seq(0.2, 5, 0.2)){
     vl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]))
     beta1 <- 1
     beta2 <- 1
     beta3 <- 0.5
                                        theta <- 1
     vl1[, 1] <- vl[1 : m]
     vl2[, 1] <- vl[ (m + 1): (m + f)]
     vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
    # vl1[, 1]<- approxfun(nvl1[, 2], nvl1[, 1])(vl1[, 2])
    # vl2[, 1]<- approxfun(nvl2[, 2], nvl2[, 1])(vl2[, 2])
    # vl3[, 1]<- approxfun(nvl3[, 2], nvl3[, 1])(vl3[, 2])
     parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
     paras <- parasall[4, ]
     crit <- 0
     broot <- c(beta1, beta2, beta3, log(theta))
     for(i in 1: 200){
         
        
        # logZ <- sapply(1:n, postlogz, theta, paras, resp)
        # dilogZ <- sum(logZ) - sum(zVec)
    
   # theta <- (bb[3*p+1])
         evl <- dfsane(c(log(vl1[, 1]), log(vl2[, 1]), log(vl3[, 1])), scoremaxnobb, method = 2, control = list(trace = TRUE), quiet = FALSE)$par#warpslv3((vl),  beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl31, parasall)## 
         vl <- exp(evl)
         vl1[, 1] <- vl[1 : m]
         vl2[, 1] <- vl[ (m + 1): (m + f)]
         vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
         #vl1[, 1]<- approxfun(nvl1[, 2], nvl1[, 1])(vl1[, 2])
         #vl2[, 1]<- approxfun(nvl2[, 2], nvl2[, 1])(vl2[, 2])
         #vl3[, 1]<- approxfun(nvl3[, 2], nvl3[, 1])(vl3[, 2])
         brootn<- multiroot(wrapscoreindv, broot, maxiter = 100, rtol = 1e-06, atol = 1e-08, ctol = 1e-08, useFortran = TRUE, positive = FALSE, jacfunc = NULL,  jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,  resp, cov, n)$root
         #broot <- multiroot(scoremaxnovl, c(1, 1, 0.5, 0), maxiter = 100, rtol = 1e-06, atol = 1e-08, ctol = 1e-08, useFortran = TRUE, positive = FALSE, jacfunc = NULL, jactype = "fullint",  verbose = FALSE, bandup = 1,banddown = 1)$root
         crit <- sum((brootn - broot)^2)
         beta1 <- broot[1]
         beta2 <- broot[2]
         beta3 <- broot[3]
         theta <- exp(broot[4])
         if(crit <= 1e-4){
             break
         }else{
             broot <- brootn
              print(crit)
         }
             
#         parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3)
        # paras <- parasall[4, ]
         #sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
         #subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
         #critn <- -margpartial(theta, subvl, paras, sumbb)
     
}
  

}
    
    # vl <- slvl31(beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, zVec)
    # vl1[, 1] <- vl[1 : m]
    # vl2[, 1] <- vl[ (m + 1): (m + f)]
    # vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
     print(theta)
     paras <- sapply(1 : n, getA, theta,   beta1, beta2, beta3,  vl1, vl2, vl3)[4, ]
     sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
     subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
     vcr <- c(vcr, -margpartial(theta, subvl, paras, sumbb))
}




testfU <- function(i, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n,  flag){
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
#A1 <- l1* resp[i, "y1"] * exp(sum(beta1 * cov[i, ]))#fA12(i, beta1, vl1, resp, cov)
#A2 <- l2 * resp[i, "y1"]* exp(sum(beta2 * cov[i, ]))#fA12(i, beta2, vl2, resp, cov)
#A3 <- l3 * (resp[i, "y2"] - resp[i, "y1"]) *exp(sum(beta3 * cov[i, ]))#fA3(i, beta3, vl3, resp, cov)
if(resp[i, "y1"]> 3){
    A1 <- 2 * (1 - exp(-3) + exp(-3) * (resp[i, "y1"] - 3)) * exp(sum(beta1 * cov[i, ]))
    A2 <- 2 * (1 - exp(-3) + exp(-3) * (resp[i, "y1"] - 3)) * exp(sum(beta2 * cov[i, ]))
    A3 <- 4 * (1 - exp(-3) + exp(-3) * (resp[i, "y2"] - 3)) * exp(sum(beta3 * cov[i, ]))-
        4 * (1 - exp(-3) + exp(-3) * (resp[i, "y1"] - 3)) * exp(sum(beta3 * cov[i, ]))
 }else{
      A1 <- 2* (1 - exp(-resp[i, "y1"])) * exp(sum(beta1 * cov[i, ]))
      A2 <- 2* (1 - exp(-resp[i, "y1"])) * exp(sum(beta2 * cov[i, ]))
      if(resp[i, "y2"]> 3){
          A3 <- 4 * (1 - exp(-3) + exp(-3) * (resp[i, "y2"] - 3)) * exp(sum(beta3 * cov[i, ])) - 4* (1 - exp(-resp[i, "y1"])) * exp(sum(beta3 * cov[i, ]))
      }else{
          A3 <- 4* (1 - exp(-resp[i, "y2"])) * exp(sum(beta3 * cov[i, ])) - 4* (1 - exp(-resp[i, "y1"])) * exp(sum(beta3 * cov[i, ]))
          
      }
      }

    A <- A1 + A2 + A3
                                        x <- cov[i, ]
B <- fB(i, theta, resp)
commd <- (1 + theta * A)
if(is.na(log(commd))){
    browser()
    print(c(A, theta))
}
                                        # print(c(commd, theta))
if(flag == 2){
                                        # browser()
    u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
    
}else{
    if(flag == 1){
        u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
        u2 <- d1 * x - B * x * theta * A1 / commd ##a p \times 1 vector
        u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
        u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
        mu <- matrix(c(u1, u2, u3, u4), nrow = 1) #(3 * p + 1)vector
        return(mu)
    }else{
        u5 <- exp(t(beta1) %*% x) #B * theta * exp(t(beta1) %*% x) / commd
        u6 <- exp(t(beta2) %*% x)#B * theta * exp(t(beta2) %*% x) / commd
        u7 <- exp(t(beta3) %*% x)#B * theta * exp(t(beta3) %*% x) / commd
        return(cbind(u5, u6, u7))
    }
}

}



testscore <- function(i, bb, resp, cov){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
    theta <- exp(bb[3*p + 1])
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- l1 * resp[i, "y1"]^(a1) * exp(sum(beta1 * cov[i, ]))#fA12(i, beta1, vl1, resp, cov)
    A2 <- l2 * resp[i, "y1"]^(a2)* exp(sum(beta2 * cov[i, ]))#fA12(i, beta2, vl2, resp, cov)
    A3 <- l3 * (resp[i, "y2"]^(a3) - resp[i, "y1"]^(a3)) *exp(sum(beta3 * cov[i, ]))#fA3(i, beta3, vl3, resp, cov)
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
    mu <- c( u2, u3, u4, u1) #(3 * p + 1)vector
    return(mu)
    


}


wraptestscore <- function(bb, resp, cov, n){
    apply(sapply(1 : n, testscore, bb, resp, cov), 1, sum)
    }
testslvtbb <- function(tbb, vl1, vl2, vl3, resp, cov, n, flag = 1){
if(flag == 1){
    theta <- exp(tbb[1]) 
    beta1 <- tbb[2 : (p + 1)]
    beta2 <- tbb[(p + 2) : (2 * p + 1)]
    beta3 <- tbb[(2 * p + 2) : (3 * p + 1)]
    lres <- lapply(1 : n, testfU,theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n,  flag = 1)
    mu <- do.call(rbind, lres)
    scu <- apply(mu, 2, sum)
}else if(flag == 2){
    theta <- (tbb[1])
    beta1 <- tbb[2 : (p + 1)]
    beta2 <- tbb[(p + 2) : (2 * p + 1)]
    beta3 <- tbb[(2 * p + 2) : (3 * p + 1)]
    lres <- sapply(1 : n, testfU,theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 2)
    scu <- sum(lres)
}
}
testwrapfun <- function(t, beta1, beta2, beta3,  vl1, vl2, vl3, resp, cov, n, flag = 2){
tbb <- c((t), beta1, beta2, beta3)
testslvtbb(tbb, vl1, vl2, vl3, resp, cov, n, flag = 2) 
}

testslvl2 <- function( theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0){
`vd1 <- resp[, "d1"]
d2 <- resp[, "d2"]
y1 <- resp[, "y1"]
y2 <- resp[, "y2"]
                                        #    subvl1 <- exp(subvl[1 : m])
                                        #   subvl2 <- exp(subvl[(m + 1) : (m + f)])
                                        #  subvl3 <- exp(subvl[(m + f + 1) : (m + f + g)])
                                        # vl1[, 1] <- subvl1 #approxfun(subvl1[, 2], subvl1[, 1])(vl1[, 2])
##vl2[, 1] <-subvl2 #approxfun(subvl2[, 2], subvl2[, 1])(vl2[, 2])
                                        # vl3[, 1] <- subvl3#approxfun(subvl3[, 2], subvl3[, 1])(vl3[, 2])
lres <- lapply(1 : n, testfU, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)
mu <- do.call(rbind, lres)
fvl <- function(j, vl, submu, flg){
    if(flg == 1){
        ix <- y1 >= vl[j, 2]
    }else 
        ix <- (y2 >= vl[j, 2]) & (y1 < vl[j, 2])
    1 /  (sum(submu[ix]) + 1e-10)
}
svl1 <- sapply(1 : m, fvl, vl1, mu[, 1], 1)
svl2 <- sapply(1 : f, fvl, vl2, mu[, 2], 1)
svl3 <- sapply(1: g, fvl, vl3, mu[, 3], 0)
return(c(svl1, svl2, svl3))

}

res <- rep(NA, nsim)
lres <- matrix(NA, nsim, 4)
for(itr  in 1 : nsim){ 
res[itr] <-  uniroot(testwrapfun, c(0.0001, 100), b1, b2, b3, vl1, vl2, vl3, lsimresp1[[itr]], lcovm[[itr]],n,  2)$root
}


for(itr in 1: nsim){
lres[itr, ]<- multiroot(wraptestscore, c( 1, 1, 0.5, 0), maxiter = 100, rtol = 1e-06, atol = 1e-08, ctol = 1e-08, useFortran = TRUE, positive = FALSE, jacfunc = NULL,  jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,  lsimresp1[[itr]], lcovm[[itr]], n)$root
}

findla<- function(la, vl){
l <- la[1]
a <- la[2]
sum((cumsum(vl[, 1]) - l * vl[, 2]^(a))^2)
}
a <- optim(c(1, 2), findla, gr = NULL, mvl1[[1]][[3]], method = "L-BFGS-B", lower = c(0.1, 0.1), upper = c(10, 10), control = list(), hessian = FALSE)
