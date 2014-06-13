library(R2Cuba)
library(MASS)
library(BB)
library(survival)
library(pracma)
hzd1 <- function(kappa1,  beta1,  t, x, g){
    exp((trans1(t)-t(beta1) %*% x-log(g))/kappa1^2) * 1/t * 1/kappa1^2
    
}

hzd2 <- function(kappa2,  beta2,  t, x, g){
    exp((trans2(t)-t(beta2) %*% x-log(g))/kappa2^2) * 1/t * 1/kappa2^2
    
}

hzd3 <- function(kappa3,  beta3,  t, x, g){
    exp((trans3(t)-t(beta3) %*% x-log(g))/kappa3^2) * 1/t * 1/kappa3^2
    
}

surv1 <- function(kappa1,  beta1, t, x, g){
   exp(- exp( (trans1(t) - log(g) - t(beta1) %*% x) / kappa1 ^2))
}

surv2 <- function(kappa2, beta2, t, x, g){
   exp(- exp( (trans2(t) - log(g) - t(beta2) %*% x) / kappa2 ^2))
}

surv3 <- function(kappa3,  beta3, t, x, g){
   exp(- exp( (trans3(t) - log(g) - t(beta3) %*% x) / kappa3 ^2))
}

ngcumhzd1 <- function(kappa1,  beta1, t, x, g){
   (- exp( (trans1(t) - log(g) - t(beta1) %*% x) / kappa1 ^2))
}

ngcumhzd2 <- function(kappa2, beta2, t, x, g){
   (- exp( (trans2(t) - log(g) - t(beta2) %*% x) / kappa2 ^2))
}

ngcumhzd3 <- function(kappa3,  beta3, t, x, g){
   (- exp( (trans3(t) - log(g) - t(beta3) %*% x) / kappa3 ^2))
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

lkhd.exp <- expression((( ((log(y1)-b1x -log(g))/kappa1^2) + log( 1/y1) + log( 1/kappa1^2)) *  d1  +  ( ((log(y1)-b2x-log(g))/kappa2^2) + log(1/y1) +log( 1/kappa2^2)) * (1-d1) +  (((log(y2)-b3x-log(g))/kappa3^2) + log( 1/y2) + log( 1/kappa3^2)) *  d1  +  (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) -(- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * d1 -  ((- exp( (log(v) - log(g) -b1x) / kappa1^2)) + (- exp( (log(v) - log(g) - b2x) / kappa2^2))) +  log(1 / (1 - vt11)) + log(1/ (1 - vt12)) ))

lkhd.exp <- expression((( ((log(y1)-b1x -log(g))/kappa1^2) + log( 1/y1) + log( 1/kappa1^2)) *  d1  +  ( ((log(y1)-b2x-log(g))/kappa2^2) + log(1/y1) +log( 1/kappa2^2)) * (1-d1) +  (((log(y2)-b3x-log(g))/kappa3^2) + log( 1/y2) + log( 1/kappa3^2)) *  d1  +  (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) -(- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * d1 +  log(1 / (1 - vt11)) + log(1/ (1 - vt12)) ))
delike <- function(y1, y2, d1, b1x, b2x, b3x,  g, kappa1, kappa2, kappa3 ){
    dd1x <- -1/kappa1^2 * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-1 / kappa1^2)
    dd2x <- -1/kappa2^2 * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-1 / kappa2^2)
    dd3x <- -1/kappa3^2 * (d1) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * (-1/kappa3^2)
    dk1 <- (-2 * (log(y1)-b1x -log(g))/kappa1^3 + (-2/kappa1)) * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-2 * (log(y1)-b1x -log(g))/kappa1^3)
    dk2 <- (-2 * (log(y1)-b2x -log(g))/kappa2^3 + (-2/kappa2)) * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-2 * (log(y1)-b2x -log(g))/kappa2^3)
    dk3 <- (-2 * (log(y2)-b3x -log(g))/kappa3^3 + (-2/kappa3)) * d1 + (- exp( (log(y2) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y2)-b3x -log(g))/kappa3^3) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y1)-b3x -log(g))/kappa3^3)
    return(cbind(dd1x, dd2x, dd3x, dk1, dk2, dk3))
}

dlike <- deriv(lkhd.exp, c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3"))
score <- function(vt,orgt,  theta, x, g, v= 1e-5){
    vt1 <- vt
    vt11 <- vt1[, 1]
    vt12 <- vt1[, 2]
    vt <- orgt
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- t(beta1)%*%x
    b2x <- t(beta2)%*%x
    b3x <- t(beta3)%*%x
    d1 <- as.numeric(vt[, 1] < vt[, 2])
    y1 <- pmin(vt[, 1], vt[, 2])
    y2 <- vt[, 2]
    
    derivlike <- delike(y1, y2, d1, b1x, b2x, b3x, g, kappa1, kappa2, kappa3)# attributes(eval(dlike))$gradient
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    score <- cbind( derivlike[, 4: ncol(derivlike)], derivlike[, 1] %*% matrix(x, ncol = p), derivlike[, 2] %*% matrix(x, ncol = p), derivlike[, 3] %*% matrix(x, ncol = p))
}

singlescore <- function(vt, orgt,  theta, x, g, v= 1e-5){
    vt1 <- vt
    vt11 <- vt1[1]
    vt12 <- vt1[2]
    vt <- orgt
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    b1x <- t(beta1)%*%x
    b2x <- t(beta2)%*%x
    b3x <- t(beta3)%*%x
    d1 <- as.numeric(vt[1] < vt[2])
    y1 <- pmin(vt[1], vt[2])
    y2 <- vt[2]
    
    derivlike <-delike(y1, y2, d1, b1x, b2x, b3x, g, kappa1, kappa2, kappa3)# attributes(eval(dlike))$gradient
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    score <- c(derivlike[4: length(derivlike)], derivlike[1] * (x), derivlike[2] * x, derivlike[3] *  x )
}
likelihood <- function(vt, orgt,  x, g,  theta, v=1e-5){
    vt1 <- vt
    vt <- orgt
    
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    
    t1 <- vt[, 1]
    t2 <- vt[, 2]
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
#    likelihood <- hzd1(kappa1,   beta1,  y1, x, g)^(d1) * hzd2(kappa2,  beta2,  y1, x, g)^(1-d1) * hzd3(kappa3,  beta3,  y2, x, g)^d1 * exp(ngcumhzd1(kappa1,   beta1,  y1, x, g) + ngcumhzd2(kappa2,  beta2,  y1, x, g) + ngcumhzd3(kappa3,  beta3,  y2, x, g) - ngcumhzd3(kappa3,  beta3,  y1, x, g)   -ngcumhzd1(kappa1,   beta1,  v, x, g) - ngcumhzd2(kappa2,  beta2,  v, x, g)) * apply(1 / (1 - vt1), 1, prod)
      likelihood <- hzd1(kappa1,   beta1,  y1, x, g)^(d1) * hzd2(kappa2,  beta2,  y1, x, g)^(1-d1) * hzd3(kappa3,  beta3,  y2, x, g)^d1 * exp(ngcumhzd1(kappa1,   beta1,  y1, x, g) + ngcumhzd2(kappa2,  beta2,  y1, x, g) + ngcumhzd3(kappa3,  beta3,  y2, x, g) - ngcumhzd3(kappa3,  beta3,  y1, x, g)   )  *  apply(exp(vt), 1, prod)
    if(sum(is.nan(likelihood)) > 0)
        browser()
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

singlelikelihood <- function(vt, orgt,  x, g,  theta, v=1e-5){
    vt1 <- vt
    vt <- orgt
    
    kappa1 <- (theta[1])
    kappa2 <- (theta[2])
    kappa3 <- (theta[3])
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    
    t1 <- vt[1]
    t2 <- vt[2]
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
    #likelihood <- hzd1(kappa1,   beta1,  y1, x, g)^(d1) * hzd2(kappa2,  beta2,  y1, x, g)^(1-d1) * hzd3(kappa3,  beta3,  y2, x, g)^d1 * exp(ngcumhzd1(kappa1,   beta1,  y1, x, g) + ngcumhzd2(kappa2,  beta2,  y1, x, g) + ngcumhzd3(kappa3,  beta3,  y2, x, g) - ngcumhzd3(kappa3,  beta3,  y1, x, g)   -ngcumhzd1(kappa1,   beta1,  v, x, g) - ngcumhzd2(kappa2,  beta2,  v, x, g))  * prod(1 / (1 - vt1))
    likelihood <- hzd1(kappa1,   beta1,  y1, x, g)^(d1) * hzd2(kappa2,  beta2,  y1, x, g)^(1-d1) * hzd3(kappa3,  beta3,  y2, x, g)^d1 * exp(ngcumhzd1(kappa1,   beta1,  y1, x, g) + ngcumhzd2(kappa2,  beta2,  y1, x, g) + ngcumhzd3(kappa3,  beta3,  y2, x, g) - ngcumhzd3(kappa3,  beta3,  y1, x, g)   )  * prod(exp(vt))
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

Amatx <- function(ij, vg, vq, theta,  x, v = 0){
    #print("A")
    i <- ij[1]
    j <- ij[2]
    dm <- function(k, vg, vt, x, theta, v){
        likelihood(vt, -log(1 - vt),  x, vg[k], theta, v)* vq[k]
    }
    A <- function(vt){
        #dnom <- (apply((sapply(1 : m, dm, vg, vt, x,   theta, v)), 1, sum) + 1e-200)
        dnom <<- (apply(sapply(1 : m, dm, vg, vt, x,   theta, v), 1, sum) )
        likelihood(vt, -log(1 - vt), x, vg[j], theta, v)* vq[j] *  likelihood(vt,-log(1 - vt),  x, vg[i], theta, v)/dnom
    }
    Aij <- my2d(A, mgrid, 1)
    mA[i, j] <- Aij
    mA <<- mA
    if(i == j){
        num <- function(k, vg,  vt, x, theta, v= v){
            score( vt,-log(1 - vt),  theta, x, vg[k], v) * vq[k]
        }
        dm <- function(k, vg, vt, x, theta, v= v){
            likelihood(vt,-log(1 - vt),  x, vg[k], theta, v)* vq[k]
        }
        b <- function(vt){
            Reduce('+', lapply(1 :m, num, vg,  vt, x,  theta, v))/(matrix(rep(dnom, q), ncol = q) ) *  matrix(rep(likelihood(vt,-log(1 - vt),  x, vg[i], theta, v), q), ncol = q)
        }
        bi <- my2d(b, mgrid,  q)#vegas (2, leth(theta) -3, b, lower = c(0.01, 0.01), upper = c(0.99,  0.99), abs.tol = 0.01)$value
        if(sum(is.nan(bi)) > 0){
            browser()
        }
        bi[is.nan(bi)] <- 0
        mb[, i] <<- bi
    }
    
    
    return(NULL)
}
bmatx <- function(i, vg, vq,  theta,  x, v=1e-5){
   # print("b")
   num <- function(k, vg,  vt, x, theta, v= v){
       score( vt, theta, x, vg[k], v) * vq[k]
   }
   dm <- function(k, vg, vt, x, theta, v= v){
        likelihood(vt, x, vg[k], theta, v)* vq[k]
    }
   b <- function(vt){
       Reduce('+', lapply(1 :m, num, vg,  vt, x,  theta, v))/(matrix(rep(apply(sapply(1 : m, dm, vg, vt, x,   theta, v), 1, sum), q), ncol = q)) *  matrix(rep(likelihood(vt, x, vg[i], theta, v), q), ncol = q)
   }
    
    #bi <-  apply(b(vtg), 2, mean)
   bi <- my2d(b, mgrid,  q)#vegas (2, leth(theta) -3, b, lower = c(0.01, 0.01), upper = c(0.99,  0.99), abs.tol = 0.01)$value
    if(sum(is.nan(bi)) > 0){
        browser()
    }
   bi[is.nan(bi)] <- 0
   mb[, i] <<- bi
   return(NULL)
}
singleAmatx <- function(ij, vg, vq, theta,  x, v = 0){
    #print("A")
    i <- ij[1]
    j <- ij[2]
    dm <- function(k, vg, vt, x, theta, v){
        singlelikelihood(vt, x, vg[k], theta, v)* vq[k]
    }
    A <- function(vt){
        singlelikelihood(vt, x, vg[j], theta, v)* vq[j] *  singlelikelihood(vt, x, vg[i], theta, v)/ sum(sapply(1 : m, dm, vg, vt, x,   theta, v))
    }
    
   
    #Aij <- area* mean(A(vtg))#mean(apply(vt, 1, A), na.rm = T)
    #mA[i, j] <<- Aij
    Aij <- vegas (2, 1, A, lower = c(0.01, 0.01), upper = c(0.99,  0.99), abs.tol = 0.01)$value
    mA[i, j] <<- Aij
    return(NULL)
}
singlebmatx <- function(i, vg, vq,  theta,  x, v=1e-5){
   # print("b")
   num <- function(k, vg,  vt, x, theta, v= v){
       singlescore( vt, theta, x, vg[k], v) * vq[k]
   }
   dm <- function(k, vg, vt, x, theta, v= v){
        singlelikelihood(vt, x, vg[k], theta, v)* vq[k]
    }
   b <- function(vt){
       Reduce('+', lapply(1 :m, num, vg,  vt, x,  theta, v))/sum(sapply(1 : m, dm, vg, vt, x,   theta, v)) *  singlelikelihood(vt, x, vg[i], theta, v)
   }
    
    #bi <- area * apply(b(vtg), 2, mean)
   bi <- vegas (2, length(theta) -3, b, lower = c(0.01, 0.01), upper = c(0.99,  0.99), abs.tol = 0.01)$value
   mb[, i] <<- bi
   return(NULL)
}
projscore <- function(vg, vq,  theta, vt, orgt,  x, a, v= 0){
    num <- function(k, vg,  vt, orgt, x, theta, v){
       (singlescore(vt, orgt, theta, x, vg[k], v) - a[, k]) *  singlelikelihood(vt, orgt,  x, vg[k], theta, v) * vq[k] 
   }
    dm <- function(k, vg, vt, orgt,  x, theta, v){
        singlelikelihood(vt, orgt, x, vg[k], theta, v)* vq[k]
    }
    apply(sapply(1 : m, num, vg,  vt, orgt,  x,  theta, v), 1, sum)/sum(sapply(1 : m, dm, vg, vt, orgt, x,   theta, v)) 
}

creata <- function(i, theta, cmptresp,  p,  mx, cmptv){
    apply(ij, 1,  Amatx, vg, vq, theta, c(1, mx[, i]), v = cmptv[i])
   # lapply(1 :m, bmatx, vg, vq, theta, mx[i,], v= cmptv[i])
    invA <- try(ginv(mA))
    if(class(invA) == "try-error"){
        browser()
    }
    a <- t(invA %*% t(mb))
    
    
}

completescore <- function(i, theta, cmptresp,orgresp,  cn, p, ma,  cmptcovm, cmptv){
    ## apply(ij, 1,  Amatx, vg, vq, theta, cmptcovm[i, ], v = cmptv[i])
    ## lapply(1 :m, bmatx, vg, vq, theta, cmptcovm[i,], v= cmptv[i])
    ## invA <- try(ginv(mA))
    ## if(class(invA) == "try-error"){
    ##     browser()
    ## }
    ## a <- t(invA %*% t(mb))
    a <- ma[[which(apply(matrix(mx, nrow = p), 1, "==", cmptcovm[i, -1, drop = F]))]]
    pjscore <-  projscore(vg, vq, theta, cmptresp[i,c("y1", "y2")], orgresp[i, c("y1", "y2")],  cmptcovm[i, ], a,  v = cmptv[i])
    if(is.nan(sum(pjscore))){
        browser()
    }
    pjscore
    
    
}
missingscore <- function(i, theta, missresp, cmptresp, mn,   p, misscovm, cmptcovm, cmptscore, cendis, missv ){
    if(missresp[i, "d1"] == 1 & missresp[i, "d2"] == 0){
        cn <- missresp[i, "y2"]
        y1 <- missresp[i, "y1"]
        x <- misscovm[i, ]
        ix <- cmptresp[, "y2"] >= cn & cmptresp[, "y1"] < cn
        if(sum(ix) > 0){
            nr <- nrow(cmptscore[ix,, drop  = F ])
            missscore <- try(apply(diag(kert(y1, cmptresp[ix, "y1"],  ht), nr, nr) %*% diag(kerx(x[-1], cmptcovm[ix, -1, drop = F], hx), nr, nr) %*%cmptscore[ix, , drop  = F ]%*% diag(cendis, nr, nr) , 2, sum) / max(sum( kert(y1, cmptresp[ix, "y1"], ht) * kerx(x[-1], cmptcovm[ix, , drop = F], hx) * cendis ),  1e-200))}
        else
            missscore <- rep(0, q)
    }else if(missresp[i, "d1"] == 0 & missresp[i, "d2"] == 0) {
        cn <- missresp[i, "y2"]
        y1 <- missresp[i, "y1"]
        x <- misscovm[i, ]
        ix <- cmptresp[, "y1"] >= cn
        if(sum(ix) > 0){
            nr <- nrow(cmptscore[ix,, drop  = F ])
            missscore <- try(apply( diag(kerx(x[-1], cmptcovm[ix, -1, drop = F ], ht), nr, nr) %*%  cmptscore[ix, , drop = F] %*% diag(cendis, nr, nr), 2, sum) / (max(sum(   kerx(x[-1], cmptcovm[ix, -1, drop =  F], hx) * cendis), 1e-200) ))
            }
        else
            missscore <- rep(0, q)
    }
    if(class(missscore) == "try-error"){
        browser()
    }
    return(missscore)
}
estm1 <- function(theta, resp, survData, covm, n, p, mv = rep(1e-5, n)){
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    cmptix <- resp[, "d2"] == 1
    covm <- matrix(covm, n, p)
    cn <- sum(cmptix)
    missix <- resp[, "d2"] == 0
    mn <- sum(missix)
    cmptresp <- resp[cmptix, ]
    cmptcovm <- covm[cmptix, , drop = F]
    cmptv <- mv[cmptix]
    missresp <- resp[missix, ]
    misscovm <- covm[missix, , drop = F]
    missv <- mv[missix]
    cmptscore <- do.call(rbind, lapply(1 : cn,completescore, theta, cmptresp, cn, p, cmptcovm, cmptv))
    #browser()
    missscore <- do.call(rbind, lapply(1 : mn, missingscore, theta, missresp, cmptresp, mn,  p, misscovm, cmptcovm, cmptscore, missv))
    #browser()
    score <- sum((apply(rbind(cmptscore, missscore), 2, sum)   )^2)
    
    
}
estm <- function(theta, resp, survData, covm, n, p, mv = rep(1e-5, n)){
    print(theta)
    theta <- c(abs(theta[1 : 3]) + 1, theta[4:q])
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    colnames(survData) <- c("y1", "d1", "y2", "d2")
    cmptix <- resp[, "d2"] == 1
    covm <- matrix(covm, n, p)
    cn <- sum(cmptix)
    missix <- resp[, "d2"] == 0
    mn <- sum(missix)
    cmptresp <- resp[cmptix, ]
    cmptcovm <- covm[cmptix, , drop = F]
    cmptv <- mv[cmptix]
    missresp <- resp[missix, ]
    misscovm <- covm[missix, , drop = F]
    missv <- mv[missix]
    ma <- lapply(1 : length(mx), creata,  theta, cmptresp,  p,  mx, cmptv)
    cmptscore <- do.call(rbind, lapply(1 : cn,completescore, theta, cmptresp, survData,  cn, p, ma,  cmptcovm, cmptv))
#    browser()
    if(mn  > 0){
        cendis <- summary(survfit(Surv(survData[, "y2"], 1 - survData[, "d2"] ))~1, times= survData[, "y2"])$surv
        missscore <- do.call(rbind, lapply(1 : mn, missingscore, theta, missresp, cmptresp, mn,  p, misscovm, cmptcovm, cmptscore, cendis[cmptix],   missv))
    #browser()
        score <- apply(rbind(cmptscore, missscore), 2, sum)
        }else{
            score <- apply((cmptscore), 2, sum)
        }
    score <- score/n
    
    
}

simuRsk <- function(i, n, p,  theta,  cen1, cen2 ,covm = NULL){
    if(is.null(covm)){
        covm <-  matrix(c(1, rbinom(1, 1, 0.5)), p, 1 )
    }
    kappa1 <- (abs(theta[1]) + 1)^2
    kappa2 <- (abs(theta[2]) + 1)^2
    kappa3 <- (abs(theta[3]) + 1)^2
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- covm
    g <- rlnorm(1, 0, 1)  
    lb1 <- exp((- t(beta1)%*%x - log(g))/kappa1)
    lb2 <- exp((- t(beta2)%*%x - log(g))/kappa2)
    lb3 <- exp((- t(beta3)%*%x - log(g))/kappa3)
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

kert <- function(t1, vt2, h){
    dnorm((t1 - vt2)/h)
}
kerx <- function(x1, vx2, h){
    vx2 <- matrix(vx2, ncol = p -1)
    x1 <- matrix(rep(x1, nrow(vx2)), ncol = p-1, byrow = T)
    h <- matrix(rep(h, nrow(vx2)), ncol = p-1, byrow = T)
    apply(dnorm((vx2 - x1)/h), 1, prod)
}



findint <- function(vv){
    set.seed(2014)
    va <- vv[1]
    vb <- vv[2]
    #vc <- vv[3]
    #vd <- vv[4]
    vt <- cbind(runif(10000, va, 1), runif(10000, vb, 1)) 
    res <- try(sum(abs(apply(score(vt, theta, covm[1], vg[1]) / (dunif(vt[,1], va, 1) * dunif(vt[, 2], vb, 1)), 2, mean) - d$integral)))
    if(is.nan(res)){
        browser()
        }
    return(res)
    
}
    
my2d <- function (f, mgrid,  nf,  ...) 
{
    
    fun <- match.fun(f)
    f <- function(vt) fun(vt, ...)
    
    mZ <- as.matrix(f(cbind(as.vector(mgrid$X), as.vector(mgrid$Y))), ncol = nf)
    temp <- function(i){
        Z <- matrix(mZ[, i], ng, ng)
        Q <- c( wx %*% Z %*% as.matrix(wy))
    }
    Q <- sapply(1 : nf, temp)
    return(Q)
}
######################################################



m = 15
theta <- c(0.3, 0.3, 0.3,  0.1, 0.5, 0.8, 0.6, 0.5,  0.5)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 500
 p <- 2
ij <- as.matrix(expand.grid(1 : m, 1 : m))
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 32
cx <- gaussLegendre(ng, 0.01, 0.99)
x <- cx$x
wx <- cx$w
cy <- gaussLegendre(ng, 0.01, 0.99)
y <- cy$x
wy <- cy$w
mgrid <- meshgrid(x, y)    
set.seed(2015)
simall <- function(itr){
    set.seed(itr + 2014)
    do.call(rbind, lapply(1:n, simuRsk, n, p, theta, 300, 400))
}
#lsurvData <- lapply(1 : 100, simall)
evalestm <- function(itr){
    survData <- lsurvData[itr]
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    ht <- n^(-1/3) * bw.nrd0(resp[resp[, "d1"] == 1, "y1"])
    hx <- n ^ (-1/3) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    vg <- quantile(survData[, 5 + p], seq(0, 1, length.out = m))
    vq <- dlnorm(vg, 0, 1) /sum(dlnorm(vg, 0, 1) )
    dfsane(theta, estm, method = 2, control = list(tol = 1.e-3, noimp = 10 ), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n))$par
}
#res <- lapply(1:10, eval)
survData <- do.call(rbind, lapply(1:n, simuRsk, n, p, theta, 300, 400))
resp <- cbind(pmin(1 - exp(-(survData[, 1])^(1/1)), 1), survData[, 2], pmin(1 - exp(-(survData[, 3])^1), 1), survData[, 4])
colnames(resp) <- c("y1", "d1", "y2", "d2")
covm <- matrix(survData[, 5 : (4+ p)], n, p)
ht <- n^(-1/3) * bw.nrd0(resp[resp[, "d1"] == 1, "y1"])
hx <- n ^ (-1/3) * apply(covm[, -1, drop = F], 2, bw.nrd0)
vg <- quantile(survData[, 5 + p],  seq(0, 1, length.out = m))
vq <- dlnorm(vg, 0, 1) /sum(dlnorm(vg, 0, 1) )

              #dfsane(theta, estm, method = 2, control = list(tol = 1.e-3, noimp = 10 ), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n))
mx <- matrix(c(0, 1), ncol = p)
#multiroot(estm, c(rep(1, 6), rep(-0.5, 3)), maxiter = 100,  rtol = 1e-6, atol = 1e-8, ctol = 1e-8,useFortran = TRUE, positive = FALSE,jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,resp,survData[, 1:4],  covm, n, p)
#theta<- c(0.08241974,  -0.06403001,   0.21495395,   0.50000000,   0.50000000,   0.50000000,  -0.44144138,  -0.50645970,  -0.85097759)

rootfind <- function(theta, otheta){
    theta <- solve(jacobian(estm, otheta, method="simple", resp, covm, n, p)) %*% estm(otheta, resp, covm, n, p) + otheta
}
    
    
delike1 <- function(y1, y2, d1, b1x, b2x, b3x,  g, kappa1, kappa2, kappa3 ){
    dd1x <- -1/kappa1^2 * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-1 / kappa1^2)
    dd2x <- -1/kappa2^2 * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-1 / kappa2^2)
    dd3x <- -1/kappa3^2 * (d1) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * (-1/kappa3^2)
    dk1 <- (-2 * (log(y1)-b1x -log(g))/kappa1^3 + (-2/kappa1)) * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-2 * (log(y1)-b1x -log(g))/kappa1^3)
    dk2 <- (-2 * (log(y1)-b2x -log(g))/kappa2^3 + (-2/kappa2)) * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-2 * (log(y1)-b2x -log(g))/kappa2^3)
    dk3 <- (-2 * (log(y2)-b3x -log(g))/kappa3^3 + (-2/kappa3)) * d1 + (- exp( (log(y2) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y2)-b3x -log(g))/kappa3^3) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y1)-b3x -log(g))/kappa3^3)
    return(cbind(dd1x, dd2x, dd3x, dk1 +  dk2+ dk3)) #
}
vsinglescore <-  function(vt, orgt, theta, x, g, v= 1e-5){
    vt1 <- vt
    vt11 <- vt1[, 1]
    vt12 <- vt1[, 2]
    vt <- orgt
    kappa1 <- (theta[1])
    kappa2 <- (theta[1])
    kappa3 <- (theta[1])
    beta1 <- theta[2 : (1 + p)]
    beta2 <- theta[(2 + p) : (1 + 2 * p)]
    beta3 <- theta[(2 + 2* p) : (1 + 3 * p)]
    b1x <- x %*% beta1
    b2x <- x %*% beta2
    b3x <- x %*% beta3
    d1 <- as.numeric(vt[, 1] < vt[, 2])
    y1 <- pmin(vt[, 1], vt[, 2])
    y2 <- vt[, 2]
    
    derivlike <- delike1(y1, y2, d1, b1x, b2x, b3x, g, kappa1, kappa2, kappa3)# attributes(eval(dlike))$gradient
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    score <- cbind( derivlike[, 4], diag(derivlike[, 1]) %*% matrix(x, ncol = p), diag(derivlike[, 2]) %*% matrix(x, ncol = p), diag(derivlike[, 3]) %*% matrix(x, ncol = p))
}
estm2 <- function(theta, resp, covm,  n, p){
    theta <- c(abs(theta[1]) + 1, theta[2 : (q - 2)])
    apply(( vsinglescore(resp[, c(1, 3)], survData[, c(1, 3)], theta, covm, survData[, 5+ p])), 2,  sum)
}

#temp<- dfsane(c(theta[1], theta[4:q]), estm2, method = 2, control = list(tol = 1.e-7, noimp = 100 ), quiet = FALSE, resp, covm, n, p)
estm(theta, resp, survData[, 1:4], covm, n, p)
