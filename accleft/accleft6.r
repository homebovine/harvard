
                                        #library(R2Cuba)
library(MASS)
library(BB)
library(survival)
#library(pracma)
library(parallel)
hzd1 <- function(kappa1,  beta1,  t, x, g){
    g * exp((trans1(t)-t(beta1) %*% x)/kappa1) * 1/t * 1/kappa1
    
}

hzd2 <- function(kappa2,  beta2,  t, x, g){
    g* exp((trans2(t)-t(beta2) %*% x)/kappa2) * 1/t * 1/kappa2
    
}

hzd3 <- function(kappa3,  beta3,  t, x, g){
    g * exp((trans3(t)-t(beta3) %*% x)/kappa3) * 1/t * 1/kappa3
    
}

surv1 <- function(kappa1,  beta1, t, x, g){
   exp(- g * exp( (trans1(t)  - t(beta1) %*% x) / kappa1 ))
}

surv2 <- function(kappa2, beta2, t, x, g){
   exp(- g * exp( (trans2(t) - t(beta2) %*% x) / kappa2 ))
}

surv3 <- function(kappa3,  beta3, t, x, g){
   exp(- g * exp( (trans3(t)  - t(beta3) %*% x) / kappa3 ))
}

ngcumhzd1 <- function(kappa1,  beta1, t, x, g){
   (- g* exp( (trans1(t)  - t(beta1) %*% x) / kappa1 ))
}

ngcumhzd2 <- function(kappa2, beta2, t, x, g){
   (- g * exp( (trans2(t)  - t(beta2) %*% x) / kappa2 ))
}

ngcumhzd3 <- function(kappa3,  beta3, t, x, g){
   (- g * exp( (trans3(t) - t(beta3) %*% x) / kappa3 ))
}

trans1 <- function(t){
    log(t)
}

trans2 <- function(t){
    log(t)
}

trans3 <- function(t){
    log(t)
}




clkhd.intg <- expression((g * exp((log(y1)-b1x)/kappa1 ) * 1/y1 * 1/kappa1 ) ^ d1* (g * exp((log(y1)-b2x)/kappa2 ) * 1/y1 * 1/kappa2 ) ^ ((1 - d1) * d2) * (g * exp((log(y2)-b3x)/kappa3 ) * 1/y2 * 1/kappa3 ) ^ (d1*d2) * exp(-  g * (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))  )) )#* dbinom(x, 1, 0.8)

cllgk.intg <- expression((log(g) +  ((log(y1)-b1x)/kappa1 ) + log( 1/y1) + log( 1/kappa1) ) *  d1 +  (log(g) + ((log(y1)-b2x)/kappa2 ) + log( 1/y1) + log(1/kappa2) ) * ((1 - d1) * d2) + (log(g) + ((log(y2)-b3x)/kappa3 ) +log(1/y2) + log(1/kappa3) ) * (d1*d2) + (-  g * (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))  )) )#* dbinom(x, 1, 0.8)

lkhd.intg <- expression((exp((log(y1)-b1x)/kappa1 ) * 1/y1 * 1/kappa1 ) ^ d1* (exp((log(y1)-b2x)/kappa2 ) * 1/y1 * 1/kappa2 ) ^ ((1 - d1) *d2)* (exp((log(y2)-b3x)/kappa3 ) * 1/y2 * 1/kappa3 ) ^ (d1*d2) * (1 + nu) ^ (d1*d2) * (1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))  )) ^ (-d1 - d2 - 1/ nu ) )#* dbinom(x, 1, 0.8)

llgk.intg <- expression((((log(y1)-b1x)/kappa1 ) +log(1/y1) + log (1/kappa1 ))  *  d1 +  (((log(y1)-b2x)/kappa2 ) +  log(1/y1) + log(1/kappa2 )) * ((1 - d1)* d2) + (((log(y2)-b3x)/kappa3 ) + log(1/y2) + log(1/kappa3 )) * (d1 * d2) +  log(1 + nu) * (d1 * d2) +  log(1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 )))  ) * (-d1 - d2 - 1/ nu ))

llgk.intg1 <- expression((((log(y1)-b1x)/kappa1 ) +log(1/y1) + log (1/kappa1 ))  *  d1 +  (((log(y1)-b2x)/kappa2 ) +  log(1/y1) + log(1/kappa2 )) * ((1 - d1) * d2) + (((log(y2)-b3x)/kappa3 ) + log(1/y2) + log(1/kappa3 )) * (d1*d2) +  log(1 + nu) * (d1 * d2) +  log(1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 )))  ) * (-d1 - d2 - 1/ nu ))

dllgk.intg <- deriv(llgk.intg,  c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3"))
dcllgk.intg <- deriv(cllgk.intg,  c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3"))
dllgk.intg1 <- deriv(llgk.intg1,  c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3", "nu"))
score <- function(vt,  theta, x, v= 1e-5){
    vt <- matrix(vt, ncol = 3)
    x <- matrix(x, ncol = p)
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- x %*% beta1#t(beta1)%*%x
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    
    
    derivlike <- attributes(eval(dllgk.intg))$gradient
    derivlike <- matrix(derivlike, ncol = 6)
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    l <- nrow(x)
    score <- cbind( derivlike[, 4: ncol(derivlike), drop = F], diag(derivlike[, 1], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 2], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 3], l, l) %*% matrix(x, ncol = p))
}

scoregvg <- function(vt, g,   theta, x, v= 1e-5){
    vt <- matrix(vt, ncol = 3)
    x <- matrix(x, ncol = p)
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- x %*% beta1#t(beta1)%*%x
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    
    
    derivlike <- attributes(eval(dcllgk.intg))$gradient
    derivlike <- matrix(derivlike, ncol = 6)
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    l <- nrow(x)
    score <- cbind( derivlike[, 4: ncol(derivlike), drop = F], diag(derivlike[, 1], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 2], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 3], l, l) %*% matrix(x, ncol = p))
}



likelihood2 <- function(vt,    x,  theta, v=1e-5){
    vt <- matrix(vt, ncol = 3)
    x <- matrix(x, ncol = p)
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- x %*% beta1
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    t1 <- vt[, 1]
    t2 <- vt[, 2]
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    likelihood <- eval(lkhd.intg)
    
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

likelihood <- function(vt,   x, g,  theta, v=1e-5){
   
    vt <- matrix(vt, ncol = 3)
    x <- matrix(x, ncol = p)
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- x %*% beta1
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    t1 <- (vt[, 1])
    t2 <- (vt[, 2])
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    likelihood <- eval(clkhd.intg)


    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

lintg  <- function(lg, ug, vt,   x,   theta, v=1e-5){
   
    vt <- matrix(vt, ncol = 3)
    
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- matrix(x, ncol = p) %*% beta1
    b2x <- matrix(x, ncol = p) %*% beta2
    b3x <- matrix(x, ncol = p) %*% beta3
    t1 <- (vt[, 1])
    t2 <- (vt[, 2])
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    a = 1/nu + d1 + d2
    b = 1/nu +(exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))
    res <-  try((pgamma(ug, a, b) - pgamma(lg, a, b)))#eval(lkhd.intg) *
    ## if(is.nan(res)){
    ##     browser()
    ##     }
    return(res)
    
}

lkgvg <- function(g,  vt,   x,   theta, v=1e-5){
   
    vt <- matrix(vt, ncol = 3)
    
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- matrix(x, ncol = p) %*% beta1
    b2x <- matrix(x, ncol = p) %*% beta2
    b3x <- matrix(x, ncol = p) %*% beta3
    t1 <- (vt[, 1])
    t2 <- (vt[, 2])
    d1 <- as.numeric(vt[, 1] < pmin(vt[, 2], vt[, 3]))
    d2 <- as.numeric(vt[, 2] < vt[, 3])
    y2 <- pmin(vt[, 2], vt[, 3])
    y1 <- pmin(vt[, 1], y2)
    res <- eval(clkhd.intg)
    return(res)
    
}


Amatx <- function(ij, vg, vq, theta,  v = 0){
    #print("A")
    i <- ij[1]
    j <- ij[2]
    
    A <- function(mvt){
        vt <- mvt[, 1:3]
        x <- mvt[, 4:ncol(mvt)]
  
        lintg(vg[j], vg[j + 1], vt, x,  theta, v ) * likelihood(vt,  x, vg[i], theta, v)/dnom
    }
    
    Aij <- my2d1(A, cbind(XY, covm1), weight, 1)
    mA[i, j] <- Aij
    mA <<- mA
    if(i == j){
        num <- function(k, vg,  vt, x, theta, v= v){
            score( vt,  theta, x, vg[k], v) * vq[k]
        }
        dm <- function(k, vg, vt, x, theta, v= v){
            likelihood(vt,  x, vg[k], theta, v)* vq[k]
        }
        b <- function(mvt){
            vt <- mvt[, 1:3]
            x <- mvt[, 4:ncol(mvt)]
            score( vt,  theta, x, v)   *  matrix(rep(likelihood(vt,  x, vg[i], theta, v), q), ncol = q)/ matrix(rep(dnom, q), ncol = q)
        }
        
        
        bi <- my2d1(b, cbind(XY, covm1), weight, q)
        
        bi[is.nan(bi)] <- 0
        mb[, i] <<- bi
    }
    
    
    return(NULL)
}

Amatx1 <- function(ij, vg, vq, theta,  v = 0){
    #print("A")
    i <- ij[1]
    j <- ij[2]

    
fi <- lkgvg(vg[i], XY, covm1, theta)
    
    A <- function(mvt){
        vt <- mvt[, 1:3]
        x <- mvt[, 4:ncol(mvt)]
        lkgvg(vg[j], vt, x, theta) * vq[j]/ dnom * fi/dnom1
    }
    
    Aij <- my2d1(A, cbind(XY, covm1), weight, 1)
    mA[i, j] <- Aij
    mA <<- mA
    if(i == j){
        
        
        b <- function(mvt){
            vt <- mvt[, 1:3]
            x <- mvt[, 4:ncol(mvt)]
            sumscore / matrix(rep(dnom, q), ncol = q) * matrix(rep(fi, q), ncol = q) /matrix(rep(dnom1, q), ncol = q)
        }
        
        
        bi <- my2d1(b, cbind(XY, covm1), weight, q)
        
        bi[is.nan(bi)] <- 0
        mb[, i] <<- bi
    }
    
    
    return(NULL)
}

projscore <- function(vg, vq,  theta, vt, orgt,  x, a, v= 0){
    num <- function(k, vg,  vt, orgt, x, theta, v){
       ( scoregvg(vt, vg[k], theta, x, v)- a[, k]) *  matrix(rep(lkgvg(vg[k], vt, x, theta, v), q), ncol = q) * vq[k]#lintg(vg[k], vg[k + 1], vt, x, theta, v)  
   }
    
    dnom <- sum(sapply(1:m, dm, vg, vt, x, theta)) + 1e-6
    apply(sapply(1 : m, num, vg,  vt, orgt,  x,  theta, v), 1, sum)/matrix(rep(dnom, q), ncol = q)
}

creata <- function( theta, cmptresp,  p,   cmptv){
    apply(ij, 1,  Amatx1, vg, vq, theta, v = cmptv[i])
   # lapply(1 :m, bmatx, vg, vq, theta, mx[i,], v= cmptv[i])
    invA <- try(ginv(mA))
    
     a <- t(invA %*% t(mb))
    

    
}

completescore <- function(i, theta, cmptresp,orgresp,  cn, p, ma,  cmptcovm, cmptv){
    
    a <- ma
    if(cmptresp[i,c("d2")] == 1){
        centime <- cmptresp[i, "y2"] + 3
    }else{
        centime <- cmptresp[i, "y2"]
    }
    pjscore <-  projscore(vg, vq, theta, c(cmptresp[i,c("y1", "y2")], centime), orgresp[i, c("y1", "y2")],  cmptcovm[i, ], a,  v = cmptv[i])
    #if(is.nan(sum(pjscore))){
      #  browser()
    #}
    pjscore
    
    
}


missingscore <- function(i, theta, missresp, cmptresp, mn,   p, misscovm, cmptcovm, cmptscore, cendis, missv ){
    missresp <- matrix(missresp, ncol = 4)
    colnames(missresp) <- c("y1", "d1", "y2", "d2")
    if(missresp[i, "d1", drop = F] == 1 & missresp[i, "d2",  drop = F] == 0){
        cn <- as.numeric(missresp[i, "y2",  drop = F])
        y1 <- missresp[i, "y1",  drop = F]
        x <- misscovm[i, -1,  drop = F]
        ix <- cmptresp[, "y2"] >= cn & cmptresp[, "y1"] < cn
        if(sum(ix) > 0){
            nr <- nrow(cmptscore[ix,, drop  = F ])
            #missscore <- try(apply(diag(as.numeric(kert(y1, cmptresp[ix, "y1"],  ht)), nr, nr) %*% diag(as.numeric(kerx(x, cmptcovm[ix,  , drop = F], hx)), nr, nr) %*% diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( kert(y1, cmptresp[ix, "y1"], ht) * kerx(x, cmptcovm[ix, , drop = F], hx) * cendis[ix] ),  1e-200))
            #missscore <- try(apply( diag(as.numeric(kerx(x, cmptcovm[ix,  , drop = F], hx)), nr, nr) %*% diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( kerx(x, cmptcovm[ix, , drop = F], hx) * cendis[ix] ),  1e-200))
            missscore <- try(apply( diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( cendis[ix] ),  1e-200))
            
        }
        else
            missscore <- rep(0, q)
    }else if(missresp[i, "d1",  drop = F] == 0 & missresp[i, "d2",  drop = F] == 0) {
        cn <- as.numeric(missresp[i, "y2",  drop = F])
        y1 <- missresp[i, "y1",  drop = F]
        x <- misscovm[i, -1,  drop = F]
        ix <- cmptresp[, "y1"] >= cn
        if(sum(ix) > 0){
            nr <- nrow(cmptscore[ix,, drop  = F ])
            #missscore <- try(apply( diag(as.numeric(kerx(x, cmptcovm[ix, , drop = F ], hx)), nr, nr) %*% diag(cendis[ix], nr, nr)%*%  cmptscore[ix, , drop = F] , 2, sum) / (max(sum(   kerx(x, cmptcovm[ix, , drop =  F], hx) * cendis[ix]), 1e-200) ))
            missscore <- try(apply( diag(cendis[ix], nr, nr)%*%  cmptscore[ix, , drop = F] , 2, sum) / (max(sum(    cendis[ix]), 1e-200) ))
            }
        else
            missscore <- rep(0, q)
    }
    #if(class(missscore) == "try-error"){
     #   browser()
    #}
    return(missscore)
}

estm <- function(theta, resp, survData, covm, n, p, mv = rep(1e-5, n)){
    
    #print(theta)
    theta <- c(abs(theta[1 : 3])+ 0.01 , theta[4:q])
    #dnom <<- apply(sapply(1:m, dm, vg, XY, covm1, theta), 1, sum) + 1e-6
    #numerator  <- lapply(1:m, num, vg, XY, covm1, theta)
    #sumscore <<- Reduce("+", numerator)
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    colnames(survData) <- c("y1", "d1", "y2", "d2")
    cmptix <- 1:n
    
    covm <- matrix(covm, n, p)
    
    cmptresp <- resp[cmptix, ]
    cmptcovm <- covm[cmptix, , drop = F]
    cmptv <- mv[cmptix]
    ma <<- creata(theta, cmptresp,  p,   cmptv)
    cmptscore <- do.call(rbind, lapply(1 : n,completescore, theta, cmptresp, survData,  cn, p, ma,  cmptcovm, cmptv))
    score <- apply(cmptscore, 2, sum)/n
    
    
}

simuRsk <- function(i, n, p,  theta,  cen1, cen2 ,covm = NULL){
    if(is.null(covm)){
        covm <-  matrix(1, p, 1)#matrix(c(1,  rnorm(1, 0, 1)), p, 1 )
    }
    kappa1 <- (abs(theta[1]) )
    kappa2 <- (abs(theta[2]) )
    kappa3 <- (abs(theta[3]) )
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- covm
    r1 <- rbinom(1, 1, 0.6)
    
    g <- r1 * rlnorm(1, 0, 0.75) + (1 - r1) * rweibull(1, 2, scale = 2)#rgamma(1, 1/nu1, scale = nu1)  # 
    lb1 <- g * exp((- t(beta1)%*%x)/kappa1)
    lb2 <- g * exp((- t(beta2)%*%x)/kappa2)
    lb3 <- g * exp((- t(beta3)%*%x)/kappa3)
    a1 <- 1/kappa1
    a2 <- 1/kappa2
    a3 <- 1/kappa3
    p1g2 <- lb2 /(lb1 + lb2)
    r <- rbinom(1,  1, p1g2)
    c <- runif(1, cen1, cen2)
    if(r == 1){
        u <- runif(1)
        t2 <-  (- log(1- u) / ( lb1 + lb2))^(1/a2)
        t1 = t2 + 3
        
    }else{
        u <- runif(1)
        t1 <-  (- log(1- u) / (lb1 + lb2))^(1/a1)
        u <- runif(1)
        t2 <- ((-log(1 - u) + (lb3 * t1 ^ a3)) / lb3) ^ (1/a3)
    }
    y2 = min(t2, c)
    y1 = min(t1, y2)
    d1 <- as.numeric(t1 < y2)
    d2 <- as.numeric(y2 < c)
    simdata <- cbind(y1, d1, y2, d2)
    colnames(simdata) <- c("y1", "d1", "y2", "d2")
    return(c(simdata, covm, g, p1g2))
}



simuRsk3 <- function( n, p, nu,  theta,  cen1, cen2 ,covm = NULL){
    t12 <- matrix(NA, n, 3)
   
    #covm <-  matrix(c( rbinom(n, 1, 0.8)), p, 1 )
    
    kappa1 <- (abs(theta[1]) )
    kappa2 <- (abs(theta[2]) )
    kappa3 <- (abs(theta[3]) )
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- matrix(covm, ncol = p)
    lb1 <-  exp(-x %*% beta1 /kappa1) #* 1/kappa1
    lb2 <-  exp(-x %*% beta2 /kappa2)
    lb3 <-  exp(-x %*% beta3 /kappa3)
    a1 <- 1/kappa1
    a2 <- 1/kappa2
    a3 <- 1/kappa3
    p1g2 <- lb2 /(lb1 + lb2)
    r <- rbinom(n,  1, p1g2)
    ix <- r == 1
    n1 <- sum(r)
    u <- runif(n1)
    t12[ix, 2] <-   ((u^(-nu ) - 1) / (nu * (lb1[ix] + lb2[ix])))^(1/a2)
    t12[ix, 1] = t12[ix, 2] + 3
    c <- runif(n, cen1, cen2)    
    
    u <- runif(n - n1)
    t12[!ix, 1] <-   ((u^(-nu ) - 1) / (nu * (lb1[!ix] + lb2[!ix])))^(1/a1)
    u <- runif(n - n1)
    t12[!ix, 2] <-  ((u^(- nu / (1 + nu)) * ( 1 + nu * lb1[!ix] * t12[!ix, 1] ^ a1 + nu * lb2[!ix] * t12[!ix, 1]^a2) - ( 1 + nu * lb1[!ix] * t12[!ix, 1] ^ a1 + nu * lb2[!ix] * t12[!ix, 1]^a2 - nu * lb3[!ix] * t12[!ix, 1] ^ a3))/ (nu * lb3[!ix]))^ (1/a3)
    t12[, 3] <- c
    
    return(t12)
}






    
my2d <- function (f, mgrid,  nf,  ...) 
{
    
    fun <- match.fun(f)
    f <- function(vt) fun(vt, ...)
    XY <- cbind(as.vector(mgrid$X), as.vector(mgrid$Y))
    ix <- which(mgrid$X < mgrid$Y, arr.ind = TRUE)
    ix1 <- mgrid$X >= mgrid$Y
    l <- sum(ix1)
    
    mZ <- as.matrix(f(XY), ncol = nf)
    tmx <- c(wx[ix[, 1]], rep(1/l, l))
    tmy <- c(wx[ix[, 2]], which(ix1, arr.ind = TRUE)[, 2])
    temp <- function(i){
        Z <- mZ[, i]
        Q <- sum(tmx *  Z * tmy)
    }
    Q <- sapply(1 : nf, temp)
    return(Q)
}

my2d1 <- function (f, XY, weight,  nf,  ...) 
{
    
    fun <- match.fun(f)
    f <- function(vt) fun(vt, ...)
    #XY <- cbind(as.vector(mgrid$X), as.vector(mgrid$Y))
    
    mZ <- as.matrix(f(XY), ncol = nf)
    Q <- t(weight) %*% mZ
    return(Q)
}
simall <- function(itr, cen1, cen2){
    set.seed(itr + 2014)
    do.call(rbind, lapply(1:n, simuRsk, n, p,     theta, cen1, cen2))
}
evalestm <- function(itr){
    survData <- lsurvData[itr]
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    ht <- n^(-1/3) * bw.nrd0(resp[resp[, "d1"] == 1, "y1"])
    hx <- n ^ (-1/3) * apply(covm[, , drop = F], 2, bw.nrd0)
    vg <- quantile(survData[, 5 + p], seq(0, 1, length.out = m))
    vq <- dlnorm(vg, 0, 1) /sum(dlnorm(vg, 0, 1) )
    dfsane(theta, estm, method = 2, control = list(tol = 1.e-3, noimp = 10, maxit = 200), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n))$par
}

######################################################


set.seed(2014)
m = 25
m1 = 100
theta <- c(0.5, 0.5, 0.5, -0.6,    -0.3,   -0.5)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 750

p <- 1
ij <- as.matrix(expand.grid(1 : m, 1 : m))
nu1 <- 0.5
nu <- 0.5
#v <- var(do.call(rbind, lsurvData)[, 6])
#u <- mean(do.call(rbind, lsurvData)[, 6])
#shape <- u^2 /v
#scale <- v/u
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 1500
up = 20
cen1 <- 300
cen2 <- 400
mx <- matrix(c(0, 1), ncol = p)
#survData <- do.call(rbind, lapply(1:n, simuRsk, n, p,  theta, 1, 3))
#survData0 <- do.call(rbind, lapply(1:n, simuRsk1, n, p,nu,  theta0, 300, 400))
#lsurvData <- lapply(1 : 100, simall, cen1, cen2)
covm1 <<- matrix(1, ng, p)#matrix(cbind(rep(1, ng), rnorm(ng, 0, 1)), ncol = p)
    
    
weight <<- rep(1/ng, ng)    
vg1 <<- qgamma(seq(0.00000001, 0.99999999999, length.out = m1 + 1), 1/nu, 1/nu)
dm <- function(k, vg, vt, x, theta, v= v){
            lkgvg(vg[k],  vt, x, theta)* vq[k]
        }
num <- function(k, vg,  vt, x, theta, v= v){
            scoregvg(vt, vg[k],   theta, x, v= 1e-5) *matrix(rep(lkgvg(vg[k], vt, x, theta), q), ncol = q) * vq[k]
        }
       
sRoot <- function(itr){
    survData <- lsurvData[[itr]]
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    
    ht <<- sum(resp[, "d1"] == 1)^(-2/15) * bw.nrd0(log(resp[resp[, "d1"] == 1, "y1"]))
    hx <<- n ^ (-2/15) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    XY <<- simuRsk3(ng, p, nu, theta, cen1, cen2, covm1)
    
    dnom1 <<- likelihood2(XY, covm1, theta)
    vg <<- seq(quantile(lsurvData[[itr]][, 5 + p], 0), quantile(lsurvData[[itr]][, 5 + p], 1), length.out = m+1)# quantile(lsurvData[[itr]][, 5 + p], seq(0.001, 0.999, length.out= m + 1))##qlnorm(seq(0.001, 0.999, length.out = m + 1), 0, 1.5)#
    vq <<- dgamma(vg, shape = shape, scale = scale)
    vq <<- vq /sum(vq)
    dnom <<- apply(sapply(1:m, dm, vg, XY, covm1, theta), 1, sum) + 1e-6
    numerator  <- lapply(1:m, num, vg, XY, covm1, theta)
    sumscore <<- Reduce("+", numerator)
    res <- try(dfsane(theta, estm, method = 3, control = list(tol = 1.e-5, noimp = 20, maxit = 300), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)))
    print(res$convergence)
    return(c(res$par, res$residual))
}
tsRoot <- function(itr) try(sRoot(itr))
#spg(theta, estm1, gr = NULL,  project = NULL, lower = c(rep(0, 3), rep(-Inf, 3 * p)), upper = rep(Inf, 3 * p), method = 3, projectArgs= NULL, control = list(ftol = 1.e-3 ), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n))$par
#res <- mclapply(1:10, tsRoot, mc.cores = 15)

#multiroot(estm, c(rep(1, 6), rep(-0.5, 3)), maxiter = 100,  rtol = 1e-6, atol = 1e-8, ctol = 1e-8,useFortran = TRUE, positive = FALSE,jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,resp,survData[, 1:4],  covm, n, p)

    
    
delike1 <- function(y1, y2, d1, b1x, b2x, b3x,  g, kappa1, kappa2, kappa3 ){
    dd1x <- -1/kappa1^2 * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-1 / kappa1^2)
    dd2x <- -1/kappa2^2 * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-1 / kappa2^2)
    dd3x <- -1/kappa3^2 * (d1) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * (-1/kappa3^2)
    dk1 <- (-2 * (log(y1)-b1x -log(g))/kappa1^3 + (-2/kappa1)) * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-2 * (log(y1)-b1x -log(g))/kappa1^3)
    dk2 <- (-2 * (log(y1)-b2x -log(g))/kappa2^3 + (-2/kappa2)) * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-2 * (log(y1)-b2x -log(g))/kappa2^3)
    dk3 <- (-2 * (log(y2)-b3x -log(g))/kappa3^3 + (-2/kappa3)) * d1 + (- exp( (log(y2) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y2)-b3x -log(g))/kappa3^3) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y1)-b3x -log(g))/kappa3^3)
    return(cbind(dd1x, dd2x, dd3x, dk1 +  dk2+ dk3)) #
}
vsinglescore <-  function(resp, orgt, theta, x, g, v= 1e-5){
    vt <- resp[, c(1, 3)]
    d2 <- resp[, 4]
    d1 <- resp[, 2]
    vt1 <- vt
    vt11 <- vt1[, 1]
    vt12 <- vt1[, 2]
    vt <- orgt
    kappa1 <- abs(theta[1])
    kappa2 <- abs(theta[2])
    kappa3 <- abs(theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    nu <- abs(theta[4 + 3 * p])
    b1x <- x %*% beta1
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    #d1 <- as.numeric(vt[, 1] < vt[, 2])
    
    y1 <- pmin(vt[, 1], vt[, 2])
    y2 <- vt[, 2]
    
    derivlike <-  attributes(eval(dllgk.intg1))$gradient
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    score <- cbind( derivlike[, 4 : 6], diag(derivlike[, 1]) %*% matrix(x, ncol = p), diag(derivlike[, 2]) %*% matrix(x, ncol = p), diag(derivlike[, 3]) %*% matrix(x, ncol = p), derivlike[, 7])
}
estm2 <- function(theta, resp, survData, covm,  n, p){
    theta <- theta
    apply(( vsinglescore(resp, survData[, c(1, 3)], theta, covm, survData[, 5+ p])), 2,  sum)
}
#lsurvData <- mclapply(1 : 1000, simall,0.5, 1.5, mc.cores = 15)
evalestm <- function(itr){
    survData <- lsurvData[[itr]]
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    ht <- n^(-2/15) * bw.nrd0(resp[resp[, "d1"] == 1, "y1"])
    hx <- n ^ (-2/15) * apply(covm[, , drop = F], 2, bw.nrd0)
    vg <- quantile(survData[, 5 + p], seq(0, 1, length.out = m))
    vq <- dlnorm(vg, 0, 1) /sum(dlnorm(vg, 0, 1) )
    dfsane(theta1, estm2, method = 2, control = list(tol = 1.e-7, noimp = 100 ), quiet = FALSE, resp, survData, covm,  n, p)$par
}
#res <- mclapply(1 : 100, evalestm , mc.cores = 15)
                                        #res1008 <- do.call(rbind, res)
#theta1 <- c(theta, 0.5)
#estm2(theta1, resp, survData[, 1:4], covm, n, p)
## res100 <- res1008
## res250 <- res2508
## res500 <- res5008
## res750 <- res7508
## mean100 <- round(apply(res100[, 1:q], 2, median), 3)
## mean250 <- round(apply(res250[, 1:q], 2, median), 3)
## mean500 <- round(apply(res500[, 1:q], 2, median), 3)
## mean750 <- round(apply(res750[, 1:q],  2, median), 3)
## #mean1000 <- apply(res1000, 2, median)
## sd100 <- round(apply(res100[, 1:q], 2, mad), 3)
## sd250 <- round(apply(res250[, 1:q], 2, mad), 3)
## sd500 <- round(apply(res500[, 1:q], 2, mad), 3)
## sd750 <- round(apply(res750[, 1:q], 2, mad), 3)


## mse100 <- round((mean100[1:q] - theta)^2 + sd100[1:q] ^2, 4)
## mse250 <- round((mean250[1:q]-theta)^2 + sd250[1:q] ^2, 4)
## mse500 <- round((mean500[1:q] - theta)^2 + sd500[1:q] ^2, 4)
## mse750 <- round((mean750[1:q] - theta)^2 + sd750[1:q] ^2, 4)

## res <- cbind(mean100, sd100, mse100, mean250, sd250, mse250, mean500, sd500, mse500)
## pres <- apply(res, 1, paste, collapse = "&")
