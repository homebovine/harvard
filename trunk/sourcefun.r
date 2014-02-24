library("survival")
library("BB")
library(numDeriv)
library(rootSolve)
library(parallel)
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


getA <- function(i, theta,  beta1, beta2, beta3, vl1, vl2, vl3, resp, cov){
    
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
   
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)
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
        parasall <- sapply(1 : n, getA, theta,   0, 0, 0, nvl1, nvl2, nvl3, resp, covmy)
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
   
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)
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
    parasall <- sapply(1 : n, getA, theta,  0, 0, 0, vl1, vl2, vl3, resp, covmy)
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
    paras <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)[4, ]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}

scoremaxnovl <- function(vbl){
    beta1 <- vbl[( 1): (p)]
    beta2 <- vbl[( p + 1): ( 2* p)]
    beta3 <- vbl[(2* p + 1): (3* p)]
    theta <- exp(vbl[3 * p + 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)
     bsc <-  wrapslvlb2(c(beta1, beta2, beta3, theta), vl1, vl2, vl3, resp, cov, n, parasall)
  #  bsc <-  wrapslvlb2(c(beta1, beta2, beta3, log(theta)), vl1, vl2, vl3, resp, cov, n)
    -c( bsc)
   
   
}

scoremaxnobb <- function(vbl, resp, cov, beta1, beta2, beta3, theta, vl1, vl2, vl3){
    vl1<- cbind(exp(vbl[1 : m]), vl1[, 2])
    vl2<- cbind(exp(vbl[(m+1) : (m + f)]), vl2[, 2])
    vl3<- cbind(exp(vbl[(m+f + 1) : (m + f + g)]), vl3[, 2])
    #vl1[, 1]<- approxfun(nvl1[, 2], nvl1[, 1])(vl1[, 2])
    #vl2[, 1]<- approxfun(nvl2[, 2], nvl2[, 1])(vl2[, 2])
    #vl3[, 1]<- approxfun(nvl3[, 2], nvl3[, 1])(vl3[, 2])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl3, parasall)
   # bsc <-  wrapslvlb2(c(beta1, beta2, beta3, log(theta)), vl1, vl2, vl3, resp, cov, n)
    -c(vls)
   
   
}
wraplikenobb <- function(vbl){
    vl1[, 1]<- vbl[1 : m]
    vl2[, 1]<- vbl[(m+1) : (m + f)]
    vl3[, 1]<- vbl[(m+f + 1) : (m + f + g)]
    paras <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)[4, ]
    sumbb <- sum(cov1%*%  beta1)  + sum(cov2%*%  beta2)+ sum(cov3 %*%  beta3)
    subvl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    crit <- -margpartial(theta, subvl, paras, sumbb)
    
}

scoreindv <- function(i, bb, resp, cov, vl1, vl2, vl3){
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

wrapscoreindv <- function(bb, resp, cov, n, vl1, vl2, vl3){
    apply(sapply(1 : n, scoreindv, bb, resp, cov, vl2, vl2, vl3), 1, sum)
    }

wrapscoreindv1 <- function(bb, theta, resp, cov, n){
    bb <- c(bb, log(theta))
    apply(sapply(1 : n, scoreindv, bb, resp, cov), 1, mean)[1:3]
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
getfromlist<- function(itr, lst, ind){
    lst[[itr]][[ind]]
}
