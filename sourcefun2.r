
library("survival")
library("BB")
library(numDeriv)
library(rootSolve)
library(parallel)
fA12 <- function(i, beta, vl, resp, cov, n, p){
    ix <- vl[, 2] <= (resp[i, "y1"] )
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(matrix(beta, p, 1)) %*% matrix(cov[i, ], p, 1))
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov, n, p){
    ix1 <-  (vl[, 2] <= resp[i, "y2"] )
    ix2 <- (vl[, 2] <=  resp[i, "y1"])
    
    A3 <- 0
    if(sum(ix1) != 0 )
        A3 <- sum(vl[ix1, 1]) - sum(vl[ix2, 1]) 
    
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
    logZ <- digamma(1/theta + resp[i, 1]+ resp[i, 2]) - log(1/theta + A)
    return(c(A1, A2, A3, A, Z, logZ))
}



postz <- function(i, theta, paras, resp){
    (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + paras[i])
}
postlogz <- function(i, theta, paras, resp){
    digamma(1/theta + resp[i, 1]+ resp[i, 2]) - log(1/theta + paras[i])
}









slvl3 <- function( b1, b2, b3, theta, vl1, vl2, vl3, resp, cov, n, parasall, p, m, f, g, lvl1, lvl2, lvl3){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    zVec <- parasall[5, ]
    A <- parasall[4, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec 
   # browser()
    fvl <- function(j, lvl, vl,  expb, flg){
        
        a <- (sum(expb[lvl[[j]]]))
        
        if(a != 0 ){
            1/vl[j, 1] - a
        }else{
            #browser()
            0
        }
        
        
    }
    
    
    svl1 <- sapply(1 : m, fvl, lvl1, vl1,  expb1, 1) 
    svl2 <- sapply(1 : f, fvl, lvl2, vl2, expb2, 2) 
    svl3 <- sapply(1: g, fvl, lvl3, vl3,  expb3, 0) 
    
    
    return(c(svl1, svl2, svl3))
    
}
slvl31 <- function( b1, b2, b3, theta, vl1, vl2, vl3, resp, cov, n, parasall, p, m, f, g, lvl1, lvl2, lvl3){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    zVec <- parasall[5, ]
    A <- parasall[4, ]
    expb1 <- exp(cov %*% b1) * zVec
    expb2 <- exp(cov %*% b2) * zVec
    expb3 <- exp(cov %*% b3) * zVec 
    
    fvl <- function(j, lvl, vl,  expb, flg){
        
        a <- (sum(expb[lvl[[j]]]))
        
        if(a != 0 ){
            1/ a
        }else{
            #browser()
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
            ix <- y1 >= (vl[j, 2] )#&((d1 == 1) )
                                        
        }else if(flg == 2){
             ix <- (y1 >= (vl[j, 2]))# & (d1 == 0)
          
        }else{
           
            ix <- (y2 >= vl[j, 2] ) & (y1 < vl[j, 2])#& (d1 == 1)
        }
}
    

warpslv3 <- function(vl, b1, b2, b3, theta,  vl1, vl2, vl3, resp, cov, n, svl, parasall, p, m, f, g, lvl1, lvl2, lvl3){
    vl1[, 1] <- (vl[1 : m])
    vl2 [, 1] <- (vl[ (m +1) : (m + f)])
    vl3 [, 1] <- (vl[(m + f + 1) : (m + f + g)])
  #  t <- (vl[ m + f + g + 1])

    svl(b1, b2, b3, theta, vl1, vl2, vl3, resp, cov, n, parasall, p, m, f, g, lvl1, lvl2, lvl3)       
}






scoremaxnobb <- function(vbl, resp, cov, beta1, beta2, beta3,  theta, vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3, parasall){
    vl1<- cbind(exp(vbl[1 : m]), vl1[, 2])
    vl2<- cbind(exp(vbl[(m+1) : (m + f)]), vl2[, 2])
    vl3<- cbind(exp(vbl[(m+f + 1) : (m + f + g)]), vl3[, 2])
    
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    vls <- warpslv3(c((vl1[, 1]), (vl2[, 1]), (vl3[, 1])), beta1, beta2, beta3, theta,   vl1, vl2, vl3, resp, cov, n, slvl3, parasall, p, m, f, g, lvl1, lvl2, lvl3)
  
    c(vls) 
   
   
}

scoremaxnovl <- function( resp, cov, beta1, beta2, beta3,  theta, vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3, parasall){
   # 
    vls <- dfsane(c(beta1, beta2, beta3), scoreindv2, method = 2, control = list(trace = FALSE), quiet = T,  theta, resp, cov, vl1, vl2, vl3, lvl1, lvl2, lvl3, m, f, g, n, p, parasall)$par
  
    c(vls) 
   
   
}






scoreindv <- function(i, bb, theta, resp, cov, vl1, vl2, vl3, n, p, parasall){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
   # theta <- exp(bb[3 * p + 1])
    
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    
    
    zVec <- parasall[5, i]
   
    A1 <- fA12(i, beta1, vl1, resp, cov, n, p)
    A2 <- fA12(i, beta2, vl2, resp, cov, n, p)
    A3 <- fA3(i, beta3, vl3, resp, cov, n, p)
    A <- A1 + A2 + A3
    x <- cov[i, ]
    B <- fB(i, theta, resp)
    commd <- (1 + theta * A)
    ## if(is.na(log(commd))){
    ##     browser()
    ##     print(c(A, theta))
    #}
    u2 <- d1 * x - B* x * theta * A1 / commd##a p \times 1 vector
    u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
    u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
   # u2 <- d1 * x - x * A1  * zVec ##a p \times 1 vector
   # u3 <- (1 - d1) * d2 * x - x * A2  * zVec ##a p \times 1 vector
   # u4 <- d1 * d2 * x - x * A3  * zVec  ##a p \times 1 vectormu <- c(u2, u3, u4) #(3 * p + 1)vector
    

    mu <- c(u2, u3, u4)
    return(mu)
    


}

scoreindv2 <- function( bb, theta, resp, cov, vl1, vl2, vl3, lvl1, lvl2, lvl3, m, f, g, n, p, parasall){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
   # theta <- exp(bb[3 * p + 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    zVec <- parasall[5, ]
    A <- parasall[4, ]
    expb1 <- exp(cov %*% beta1) * zVec
    expb2 <- exp(cov %*% beta2) * zVec
    expb3 <- exp(cov %*% beta3) * zVec
    expb1n <- diag(as.vector(exp(cov %*% beta1) * zVec)) %*% cov
    expb2n <- diag(as.vector(exp(cov %*% beta2) * zVec)) %*% cov
    expb3n <- diag(as.vector(exp(cov %*% beta3) * zVec)) %*% cov
    
    fvl <- function(j, lvl, vl,  expb, expbn, flg){
        
       
        b <- (sum(expb[lvl[[j]]]))
         a <- apply((expbn[lvl[[j]], , drop = FALSE]), 2, sum)  /b
        if(b != 0 ){
            a
        }else{
            #browser()
            0
        }
        
        
    }
    
    
    svl1 <- matrix(sapply(1 : m, fvl, lvl1, vl1,  expb1, expb1n,  1), m, p, byrow = T) 
    svl2 <- matrix(sapply(1 : f, fvl, lvl2, vl2, expb2, expb2n, 2), f, p, byrow = T) 
    svl3 <- matrix(sapply(1: g, fvl, lvl3, vl3,  expb3, expb3n, 0), g, p, byrow = T) 
    x <- matrix(cov, n, p)
   
    u2 <- t(d1) %*% x - apply(svl1, 2, sum)#
    u3 <- t((1 - d1) * d2) %*%  x - apply(svl2, 2, sum)
    u4 <- t(d1 * d2) %*% x - apply(svl3, 2, sum)#
    

    mu <- c(u2, u3, u4)
    return(mu)
    


}
scoretheta <- function(i, theta,  beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p){
    
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov, n, p)
    A2 <- fA12(i, beta2, vl2, resp, cov, n, p)
    A3 <- fA3(i, beta3, vl3, resp, cov, n, p)
    A <- A1 + A2 + A3
    x <- cov[i, ]
    B <- fB(i, theta, resp)
 
    u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
    }
wrapstha <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p){
    sum(sapply(1 :n, scoretheta, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p))
}

scorepara <- function(i, bb, theta, resp, cov, vl1, vl2, vl3, n, p){
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1) : (2 * p)]
    beta3 <- bb[(2*p + 1) : (3*p)]
    
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
        #browser()
        print(c(A, theta))
    }
    u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
    u2 <- d1 * x - B * x * theta * A1 / commd ##a p \times 1 vector
    u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
    u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
    mu <- c(u2, u3, u4, u1) #(3 * p + 1)vector
    return(mu)
    


}

wrapscoreindv <- function(bb,  theta, resp, cov, n, vl1, vl2, vl3,  p, parasall){
    apply(sapply(1 : n, scoreindv, bb, theta,  resp, cov, vl2, vl2, vl3, n, p, parasall), 1, sum)[1 : (3 * p)]
    }
wrapscoreindv1 <- function(bb, resp, cov, n, vl1, vl2, vl3,  p){
    apply(sapply(1 : n, scorepara, bb, resp, cov, vl2, vl2, vl3, n, p), 1, sum)
    }



    
margpartial <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, cov1, cov2, cov3 ){
    vl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    paras <- parasall[4, ]
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]  
    sum(resp[, 1] * resp[, 2]) * log(theta + 1)  - sum(B * log(1 + theta* paras))+  sumbb + sum(log(vl))
 #   browser()
   # eta <- 1/theta
    #log(prod(eta^eta /gamma(eta) * gamma(eta + resp[, 1] + resp[, 2])/ (eta + paras ) ^(eta + resp[, 1] + resp[, 2])))  + sum(log(vl)) + (sumbb)
}


margpartial1 <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, cov1, cov2, cov3 ){
    vl <- c(vl1[, 1], vl2[, 1], vl3[, 1])
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    paras <- parasall[4, ]
    zVec <- parasall[5, ]
    logZ <- parasall[6, ]
    sumbb <- sum(matrix(cov1, ncol = p)%*%  matrix(beta1))  + sum(matrix(cov2, ncol = p)%*%  matrix(beta2))+ sum(matrix(cov3, ncol = p)%*%  matrix(beta3))
    B <- 1/theta + resp[, 1] + resp[, 2]  
  sum(log(dgamma(zVec, 1/theta, 1/theta))) +  sumbb + sum(log(vl)) + sum((resp[, 1] + resp[, 2]) * log(zVec)) - sum(zVec * paras)
}



getfromlist<- function(itr, lst, ind){
    lst[[itr]][[ind]]
}
getvlist <- function(i, estres, ix, ix1){
    estres[[ix]][[2]][[ix1]][, 1]
}

comscore1 <- function(vlbb, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall){
    vl <- (vlbb[1 : (m + f + g)])
    vl1[, 1] <- (vlbb[1 :m ])
    vl2[, 1] <- (vlbb[(m + 1) : (m + f) ])
    vl3[, 1] <- (vlbb[(m + f + 1) : (m + f + g) ])
    beta1 <- vlbb[( m + f + g+ 1) : (m + f + g + p) ]
    beta2 <- vlbb[(m + f + g + p + 1) : (m + f + g + 2* p) ]
    beta3 <- vlbb[(m + f + g + 2 * p + 1) : (m + f + g + 3 * p)]
    theta <- (vlbb[m + f + g + 3 * p + 1])
    nobb <- scoremaxnobb(log(vl), resp, cov, beta1, beta2, beta3, theta, vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3, parasall) 
    novl <- wrapscoreindv(c(beta1, beta2, beta3), theta, resp, cov,  n, vl1, vl2, vl3,  p, parasall)
    sth <- wrapstha(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    c(sth, novl,  nobb) 
    
}


comscore3 <- function( theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3){
 
   bb <- scoremaxnovl(resp, cov, beta1, beta2, beta3, theta, vl1, vl2, vl3, n, p, m, f, g, lvl1, lvl2, lvl3, parasall)
   beta1 <- bb[1 : p]
   beta2 <- bb[(p + 1) : (2 * p)]
   beta3 <- bb[(2 * p + 1) : (3 * p)]
   parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
   vl <- warpslv3(c(vl1[, 1], vl2[, 1], vl3[, 1]), beta1, beta2, beta3, theta, vl1, vl2, vl3, resp, cov, n, slvl31, parasall, p, m, f, g, lvl1, lvl2, lvl3)
   vl1 <- cbind(vl[1 : m], vl1[, 2])
   vl2 <- cbind(vl[ (m + 1): (m + f)], vl2[, 2])
   vl3 <- cbind(vl[ (m + f + 1): (m + f+ g)], vl3[, 2])
   crit <- margpartial(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, cov1, cov2, cov3)
  
   (list(beta1, beta2, beta3, vl1, vl2, vl3, crit)) 
    
}





objfun2 <- function(theta, rtime, tol, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3){
    mglk <- 0
    for(i in 1 : rtime){
        resbb <- comscore3(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3)
        crit <- mglk - resbb[[7]]#c(resbb[[1]] - beta1,resbb[[2]] - beta2, resbb[[3]] - beta3)
        
        beta1 <- resbb[[1]]
        beta2 <- resbb[[2]]
        beta3 <- resbb[[3]]
        vl1 <- resbb[[4]]
        vl2 <- resbb[[5]]
        vl3 <- resbb[[6]]
        mglk <- resbb[[7]]
        if(max(abs(crit)) <= tol){
            break
            }
        }
    g.beta1 <<- beta1
    g.beta2 <<- beta2
    g.beta3 <<- beta3
    g.vl1 <<- vl1
    g.vl2 <<- vl2
    g.vl3 <<- vl3
    
    partlike <- -margpartial(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, cov1, cov2, cov3)

        
    }

scoreobj2 <- function(theta, rtime, tol, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, ini, verbose){
    mglk <- 0
    for(i in 1 : rtime){
        resbb <- comscore3(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3)
        crit <- (mglk - resbb[[7]])#c(resbb[[1]] - beta1,resbb[[2]] - beta2, resbb[[3]] - beta3)
        
        
        beta1 <- resbb[[1]]
        beta2 <- resbb[[2]]
        beta3 <- resbb[[3]]
        vl1 <- resbb[[4]]
        vl2 <- resbb[[5]]
        vl3 <- resbb[[6]]
        mglk <- resbb[[7]]
        if(verbose == 1){
            print(paste("iter=", i, "reduction=", crit), sep = "")
        }else if(verbose == 2){
            print(paste("iter=", i, "reduction=", crit), sep = "")
            print(c(theta, beta1, beta2, beta3))
            }
        if(mean(abs(crit)) <= tol){
            break
        }
        }
    if(ini == 1){
        partlike <- -mglk#-margpartial(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, cov1, cov2, cov3)## wrapstha(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)^2#mglk###  ##
    }else if(ini == 2){
        #browser()
        res <- c(beta1, beta2, beta3,  theta, i, mglk)
        attributes(res) <- list("vl1" = vl1, "vl2" = vl2, "vl3" = vl3)
        return(res)
    }else{
        return(list(vl1, vl2, vl3, beta1, beta2, beta3, i, mglk))
    }

        
}

simeval <- function(sitr){
   res <-  try(estreal(rep(0, 3 * p), c(0.2, 1.8), c(0.1, 2),  lsimresp1[[sitr]], lcovm[[sitr]], FALSE, 100, 1e-5)[[1]])
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
    c <- runif(1, cen1, cen2) #censoring time
    c(t1, t2, c)
}

estreal <- function(bb, vtheta, altheta,  resp, cov, hessian, rtime, tol, verbose, factr){
    p <- ncol(cov)
    n <- nrow(resp)
    resp[, 3] <- round(resp[, 3], 8)
    resp[, 4] <- round(resp[, 4], 8)
    colnames(resp) <- c("d1", "d2", "y1", "y2")
    beta1 <- bb[1 : p]
    beta2 <- bb[(p + 1): (2 * p)]
    beta3 <- bb[(2 * p + 1): (3 * p)]
    nth <- length(vtheta)
    theta <- mean(vtheta)
    d1 <- resp[, 1]
    d2 <- resp[, 2]
    y1 <- resp[, 3]
    y2 <- resp[, 4]
    ind1 <- which(d1 == 1 )
    ind2 <- which((d1 == 0) & (d2 == 1))
    ind3 <- which((d1 ==1) & (d2 == 1))
    subgix1 <- (d1 == 0)
    respsub1 <- resp[subgix1, ]
    respsub2 <- resp[ind1, ]
    surv1 <- (coxph(Surv(resp[, "y1"], resp[, "d1"]) ~ ., as.data.frame(cov)))
    surv2 <- (coxph(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~., as.data.frame(cov[subgix1, ] )))
    surv3 <- (coxph(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~., as.data.frame(cov[ind1, ])))
    bz10 <- basehaz(surv1)
    bz20 <- basehaz(surv2)
    bz30 <- basehaz(surv3)
    ixbz1 <- which(c(bz10[1, 1], diff(bz10[, 1])) != 0) 
    ixbz2 <- which(c(bz20[1, 1], diff(bz20[, 1])) != 0)
    ixbz3 <- which(c(bz30[1, 1], diff(bz30[, 1])) != 0) 
    bz1 <-  bz10[ixbz1, ]
    bz2 <-  bz20[ixbz2, ]
    bz3 <-  bz30[ixbz3, ]
    bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
    bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
    bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))
                                        #browser()
    vl10 <- bz1
    vl20 <- bz2
    vl30 <- bz3
    m <- nrow(vl10)
    f <- nrow(vl20)
    g <- nrow(vl30)
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
    vl3 <- vl30
    cov1 <- cov[ind1, ]
    cov2 <- cov[ind2, ]
    cov3 <- cov[ind3, ]
    lvl1 <- lapply(1 : m, funcvl, vl1, 1, resp)
    lvl2 <- lapply(1 : f, funcvl, vl2, 2, resp)
    lvl3 <- lapply(1 : g, funcvl, vl3, 3, resp)
    vl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]))
    
    vl1[, 1] <- vl[1 : m]
    vl2[, 1] <- vl[ (m + 1): (m + f)]
    vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
    broot <- bb
    broot0 <- c(rep(0, 3 * p), -0.5)
    
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
    temp <- optim((vtheta[1] + vtheta[2])/2, scoreobj2, gr = NULL,  rtime, tol, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 1, verbose, method = "L-BFGS-B", lower = vtheta[1], upper = vtheta[2], control = list(fnscale = 1, factr = factr), hessian = FALSE)#try(uniroot(scoreobj2, vtheta, rtime, tol,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 1, verbose))
     if(temp$convergence != 0){
         temp <- try(optim((altheta[1] + altheta[2])/2, scoreobj2, gr = NULL,  rtime, tol, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 1, verbose, method = "L-BFGS-B", lower = altheta[1], upper = altheta[2], control = list(fnscale = 1, factr = factr), hessian = FALSE))#try(uniroot(scoreobj2, altheta, rtime, tol,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 1, verbose))
         if(temp$convergence != 0){
             return(list(rep(NA, 3 * p + 1), temp))
             }
         }
     theta  <- temp$par
     res <- scoreobj2(theta, rtime, tol, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 0, verbose)
     vl1 <- res[[1]]
     vl2 <- res[[2]]
     vl3 <- res[[3]]
     beta1 <- res[[4]]
     beta2 <- res[[5]]
     beta3 <- res[[6]]
     miter <- res[[7]]
     #browser()
     if(hessian){

         hm <- try(jacobian(comscore1, c(c(vl1[, 1], vl2[, 1], vl3[, 1]), beta1, beta2, beta3, (theta)), method = "simple", method.args= list(eps = 1e-8), vl1, vl2, vl3, resp, cov, n, p, m, f, g, lvl1, lvl2, lvl3, parasall))
         if(class(hm) == "try-error"){
             print(hm)
         }
       
         mres <- c(beta1, beta2, beta3, theta)
         mvl <- list(vl1, vl2, vl3)
         res <- list( thbb = c(theta, beta1, beta2, beta3), vv1 = vl1, vv2 =vl2, vv3 =vl3, hessian = hm, miter = miter)
         class(res) <- "FrqID"
         return(res)
     }else{
         res <- list( thbb = c(theta, beta1, beta2, beta3), vv1 = vl1, vv2 =vl2, vv3 =vl3, miter = miter)
         class(res) <- "FrqID"
         return(res)
     }
 
 }


iniestreal <- function(vtheta, bb,  resp, cov, hessian, rtime, tol, verbose){
     p <- ncol(cov)
     n <- nrow(resp)
     resp[, 3] <- round(resp[, 3], 8)
     resp[, 4] <- round(resp[, 4], 8)
     colnames(resp) <- c("d1", "d2", "y1", "y2")
     beta1 <- bb[1 : p]
     beta2 <- bb[(p + 1): (2 * p)]
     beta3 <- bb[(2 * p + 1): (3 * p)]
     nth <- length(vtheta)
     theta <- (vtheta)
     d1 <- resp[, 1]
     d2 <- resp[, 2]
     y1 <- resp[, 3]
     y2 <- resp[, 4]
     ind1 <- which(d1 == 1 )
     ind2 <- which((d1 == 0) & (d2 == 1))
     ind3 <- which((d1 ==1) & (d2 == 1))
     subgix1 <- (d1 == 0)
     respsub1 <- resp[subgix1, ]
     respsub2 <- resp[ind1, ]
     surv1 <- (coxph(Surv(resp[, "y1"], resp[, "d1"]) ~ ., as.data.frame(cov)))
     surv2 <- (coxph(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~., as.data.frame(cov[subgix1, ] )))
     surv3 <- (coxph(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~., as.data.frame(cov[ind1, ])))
     bz10 <- basehaz(surv1)
     bz20 <- basehaz(surv2)
     bz30 <- basehaz(surv3)
     ixbz1 <- which(c(bz10[1, 1], diff(bz10[, 1])) != 0) 
     ixbz2 <- which(c(bz20[1, 1], diff(bz20[, 1])) != 0)
     ixbz3 <- which(c(bz30[1, 1], diff(bz30[, 1])) != 0) 
     bz1 <-  bz10[ixbz1, ]
     bz2 <-  bz20[ixbz2, ]
     bz3 <-  bz30[ixbz3, ]
     bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
     bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
     bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))
     #browser()
     vl10 <- bz1
     vl20 <- bz2
     vl30 <- bz3
     m <- nrow(vl10)
     f <- nrow(vl20)
     g <- nrow(vl30)
     n <- nrow(resp)
     vl1 <- vl10
     vl2 <- vl20
     vl3 <- vl30
     cov1 <- cov[ind1, ]
     cov2 <- cov[ind2, ]
     cov3 <- cov[ind3, ]
     lvl1 <- lapply(1 : m, funcvl, vl1, 1, resp)
     lvl2 <- lapply(1 : f, funcvl, vl2, 2, resp)
     lvl3 <- lapply(1 : g, funcvl, vl3, 3, resp)
     vl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]))
 
     vl1[, 1] <- vl[1 : m]
     vl2[, 1] <- vl[ (m + 1): (m + f)]
     vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
     broot <- bb
     #broot0 <- c(rep(0, 3 * p), -0.5)
     
     parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, p)
     temp <- scoreobj2(theta,   rtime, tol,   beta1, beta2, beta3, vl1, vl2, vl3, resp, cov,  n, p, m, f, g, lvl1, lvl2, lvl3, parasall, cov1, cov2, cov3, 2, verbose)
     return(temp) 
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
    return(cbind(simresp1, covm))
}
FrqID <- function(survData, startValues,  stheta, wtheta, hessian = F,  miter = 100, tol = 1e-4, factr = NULL, initial = F, step = 0.01, ncores = detectCores(),  verbose){
    if(is.null(factr)){
        factr <- tol
    }
    y1 <- pmin(survData[, 1], survData[, 3])
    y2 <- survData[, 3]
    d1 <- survData[, 2]
    d2 <- survData[, 4]
    resp <- cbind(d1, d2, y1, y2)
    np <- ncol(survData)
    covmy <- matrix(survData[, (5 : np)], ncol =  np - 4)
    if(!initial){
        res <- estreal(startValues, stheta, wtheta, resp, covmy,  hessian, miter, tol, verbose, factr)
        class(res) <- "FrqID"
    }else{
        stheta <- seq(stheta[1], stheta[2], step)
        res1 <- mclapply(stheta, iniestreal,  startValues,  resp, covmy,  hessian, miter, tol, verbose, mc.cores = ncores, mc.allow.recursive = FALSE)

        mres <- do.call(rbind, res1)
        nc <- ncol(mres)
        ix <- which(mres[, nc] == max(mres[, nc]))
        res <- res1[[ix]]
        attributes(res) <- list("alres" = res1)
        class(res) <- "iniFrqID"
    }
    return(res)
}
plot.iniFrqID <- function (object){
    plot(object[, 2] ~ object[, 1], type = "l", xlab = "theta values", ylab = "score functions")
}
simu <- function(itr){
    simCpRsk(5000, p = 1, theta = 1, lambda1 = 1, lambda2 = 1, lambda3 = 0.5, kappa = 2.5, beta1 = 0.5, beta2 = 0.1, beta3 = 0.3, covm = NULL, 2.5, 3)
}
analysis <- function(itr){
    survData <- lsurvData[[itr]]
    FrqID( cbind(survData), rep(0, 3), c(0.70, 1.3), c(3.51, 3.9), tol = 1e-6, factr  = NULL, initial = T, step = 0.02, ncores = 10,  verbose =1)
}
evaltheta <- function(theta){
    FrqID( cbind(lsurvData[[1]]), rep(0, 3), c(theta, theta), c(3.51, 3.9), tol = 1e-6, factr = NULL, initial = T, step = 0.02,  verbose =2)
}

lsurvData <- lapply(1 : nsim, simu)#simCpRsk(250, p = 1, theta = 1, lambda1 = 1, lambda2 = 1, lambda3 = 1, kappa = 2, beta1 = 0.5, beta2 = 0.2, beta3 = 0.3, covm = NULL, 2.5, 3)
nrealtemp <- mclapply(seq(0.89, 0.95, 0.005), evaltheta, mc.cores = 10)
lrealtemp <- mclapply(1:nsim, analysis, mc.cores = 10)
realtemp <- FrqID( cbind(lsurvData[[1]]), rep(0, 3), c(0.87, 0.89), c(3.51, 3.9), tol = 1e-6, factr = NULL, initial = T, step = 0.005,  ncores = 10,   verbose =2)
pdf(file = "temp1.pdf")
plot(realtemp[, 2] ~realtemp[, 1], type= "l")
dev.off()

