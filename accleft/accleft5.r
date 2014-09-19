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




clkhd.intg <- expression((g * exp((log(y1)-b1x)/kappa1 ) * 1/y1 * 1/kappa1 ) ^ d1* (g * exp((log(y1)-b2x)/kappa2 ) * 1/y1 * 1/kappa2 ) ^ (1 - d1) * (g * exp((log(y2)-b3x)/kappa3 ) * 1/y2 * 1/kappa3 ) ^ d1 * exp(-  g * (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))  )) )

lkhd.intg <- expression((exp((log(y1)-b1x)/kappa1 ) * 1/y1 * 1/kappa1 ) ^ d1* (exp((log(y1)-b2x)/kappa2 ) * 1/y1 * 1/kappa2 ) ^ (1 - d1) * (exp((log(y2)-b3x)/kappa3 ) * 1/y2 * 1/kappa3 ) ^ d1 * (1 + nu) ^ d1 * (1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))  )) ^ (-d1 - 1 - 1/ nu ))

llgk.intg <- expression((((log(y1)-b1x)/kappa1 ) +log(1/y1) + log (1/kappa1 ))  *  d1 +  (((log(y1)-b2x)/kappa2 ) +  log(1/y1) + log(1/kappa2 )) * (1 - d1) + (((log(y2)-b3x)/kappa3 ) + log(1/y2) + log(1/kappa3 )) * d1 +  log(1 + nu) * d1 +  log(1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 )))  ) * (-d1 - 1 - 1/ nu ))

llgk.intg1 <- expression((((log(y1)-b1x)/kappa1 ) +log(1/y1) + log (1/kappa1 ))  *  d1 +  (((log(y1)-b2x)/kappa2 ) +  log(1/y1) + log(1/kappa2 )) * ((1 - d1) * d2) + (((log(y2)-b3x)/kappa3 ) + log(1/y2) + log(1/kappa3 )) * (d1*d2) +  log(1 + nu) * (d1 * d2) +  log(1 + nu *  (exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + d1 * (exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 )))  ) * (-d1 - d2 - 1/ nu ))

dllgk.intg <- deriv(llgk.intg,  c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3"))
dllgk.intg1 <- deriv(llgk.intg1,  c("b1x", "b2x", "b3x", "kappa1", "kappa2", "kappa3", "nu"))
score <- function(vt,  theta, x, v= 1e-5){
    vt <- matrix(vt, ncol = 2)
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
    d1 <- as.numeric(vt[, 1] < vt[, 2])
    y1 <- pmin(vt[, 1], vt[, 2])
    y2 <- vt[, 2]
    
    derivlike <- attributes(eval(dllgk.intg))$gradient
    derivlike <- matrix(derivlike, ncol = 6)
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    l <- nrow(x)
    score <- cbind( derivlike[, 4: ncol(derivlike), drop = F], diag(derivlike[, 1], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 2], l, l) %*% matrix(x, ncol = p), diag(derivlike[, 3], l, l) %*% matrix(x, ncol = p))
}

singlescore <- function(vt,   theta, x,  v= 1e-5){
    vt <- vt
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
    
    derivlike <- attributes(eval(dllgk.intg))$gradient
    derivlike[is.nan(derivlike)] <- 0
    #browser()
    score <- c(derivlike[4: length(derivlike)], derivlike[1] * (x), derivlike[2] * x, derivlike[3] *  x )
}
likelihood1 <- function(vt,    x,  theta, v=1e-5){
    vt <- matrix(vt, ncol = 2)
    x <- matrix(x, nrow = p)
    x1 <- x[, -1]
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
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
    likelihood <- eval(lkhd.intg)
    #if(sum(is.nan(likelihood)) > 0)
       # browser()
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}
likelihood2 <- function(vt,    x,  theta, v=1e-5){
    vt <- matrix(vt, ncol = 2)
    x <- matrix(x, ncol = p)
    x1 <- x[, -1]
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
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
    likelihood <- eval(lkhd.intg)
    #if(sum(is.nan(likelihood)) > 0)
    #    browser()
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

likelihood <- function(vt,   x, g,  theta, v=1e-5){
   
    vt <- matrix(vt, ncol = 2)
    x <- matrix(x, ncol = p)
    x1 <- x[, -1]
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
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
    likelihood <- eval(clkhd.intg)

    #if(sum(is.nan(likelihood)) > 0)
     #   browser()
    likelihood[is.nan(likelihood)] <- 0
    return(likelihood)
}

lintg  <- function(lg, ug, vt,   x,   theta, v=1e-5){
   
    vt <- matrix(vt, ncol = 2)
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
    d1 <-  as.numeric(t1 < t2)
    y2 <-  t2
    y1 <- pmin(t1, t2)
    a = 1/nu + d1 + 1
    b = 1/nu +(exp((log(y1)-b1x)/kappa1 ) +  exp((log(y1)-b2x)/kappa2 ) + exp((log(y2)-b3x)/kappa3 ) - exp((log(y1)-b3x)/kappa3 ))
     (pgamma(ug, a, b) - pgamma(lg, a, b))#eval(lkhd.intg) *
    
}
singlelikelihood <- function(vt,  x, g,  theta, v=1e-5){
    vt <- vt
    
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
    
    likelihood <- hzd1(kappa1,   beta1,  y1, x, g)^(d1) * hzd2(kappa2,  beta2,  y1, x, g)^(1-d1) * hzd3(kappa3,  beta3,  y2, x, g)^d1 * exp(ngcumhzd1(kappa1,   beta1,  y1, x, g) + ngcumhzd2(kappa2,  beta2,  y1, x, g) + ngcumhzd3(kappa3,  beta3,  y2, x, g) - ngcumhzd3(kappa3,  beta3,  y1, x, g)   ) 
    return(likelihood)
}

Amatx <- function(ij, vg, vq, theta, x,  v = 0){
    #print("A")
    i <- ij[1]
    j <- ij[2]
    dm <- function(k, vg, vt, x, theta, v){
        likelihood(vt, -log(1 - vt),  x, vg[k], theta, v)
    }
    #XY <<- do.call(rbind, lapply(1 :ng, simuRsk2, ng,  p,  theta,  300, 400 ,x, vg[i]))
    covm1 <- cbind(rep(1, ng), rep(x, ng))
    A <- function(mvt){
        vt <- mvt[, 1:2]
        x <- mvt[, 3:ncol(mvt)]
        #dnom <- (apply((sapply(1 : m, dm, vg, vt, x,   theta, v)), 1, sum) + 1e-200)
        #dnom <<- likelihood1(vt,   x, theta0, v)
#        likelihood(vt,  x, vg[j],  theta, v)* vq1[j] *  likelihood(vt,  x, vg[i], theta, v)/dnom
        lintg(vg[j], vg[j + 1], vt, x,  theta, v ) * likelihood(vt,  x, vg[i], theta, v)/dnom
    }
    
    Aij <- my2d1(A, cbind(XY, covm1), weight, 1)#my2d(A, mgrid, 1)# + integrate(A1, min(cy$x), max(cy$x))$value
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
            vt <- mvt[, 1:2]
            x <- mvt[, 3:ncol(mvt)]
            score( vt,  theta, x, v)   *  matrix(rep(likelihood(vt,  x, vg[i], theta, v), q), ncol = q)/ matrix(rep(dnom, q), ncol = q)
        }
        
        
        bi <- my2d1(b, cbind(XY, covm1), weight, q)#my2d(b, mgrid,  q)# + adaptIntegrate(b1, min(cy$x), max(cy$x), fDim = p)$integral
        #if(sum(is.nan(bi)) > 0){
         #   browser()
        #}
        bi[is.nan(bi)] <- 0
        mb[, i] <<- bi
    }
    
    
    return(NULL)
}

projscore <- function(vg, vq,  theta, vt, orgt,  x, a, v= 0){
    num <- function(k, vg,  vt, orgt, x, theta, v){
       ( a[, k]) *  lintg(vg[k], vg[k + 1], vt, x, theta, v)  
   }
    dm <- function(k, vg, vt, orgt,  x, theta, v){
        likelihood(vt,  x, vg[k], theta, v)* vq[k]
    }
    score(vt,  theta, x) - apply(sapply(1 : m, num, vg,  vt, orgt,  x,  theta, v), 1, sum)#/likelihood1(vt,  x, theta, v)
}

creata <- function( theta, cmptresp,  p, x,   cmptv){
    apply(ij, 1,  Amatx, vg, vq, theta,  x, v = cmptv[i])
   # lapply(1 :m, bmatx, vg, vq, theta, mx[i,], v= cmptv[i])
    invA <- try(ginv(mA))
    #if(class(invA) == "try-error"){
      #  browser()
    #}
     a <- t(invA %*% t(mb))
    ## myapprox <- function(rnum){
    ##     approx(vg[1:m], a[rnum, ], vg1, rule = 2)$y
    ## }
#   do.call(rbind, lapply(1:q, myapprox))
    
}

completescore <- function(i, theta, cmptresp,orgresp,  cn, p, ma,  cmptcovm, cmptv){
    ## apply(ij, 1,  Amatx, vg, vq, theta, cmptcovm[i, ], v = cmptv[i])
    ## lapply(1 :m, bmatx, vg, vq, theta, cmptcovm[i,], v= cmptv[i])
    ## invA <- try(ginv(mA))
    ## if(class(invA) == "try-error"){
    ##     browser()
    ## }
    ## a <- t(invA %*% t(mb))
    a <- ma
    pjscore <-  projscore(vg, vq, theta, cmptresp[i,c("y1", "y2")], orgresp[i, c("y1", "y2")],  cmptcovm[i, ], a[[cmptcovm[i, 2] + 1]],  v = cmptv[i])
    #if(is.nan(sum(pjscore))){
      #  browser()
    #}
    pjscore
    
    
}

completescore1 <- function(theta, cmptresp,orgresp,  cn, p,  cmptcovm, cmptv){
    ## apply(ij, 1,  Amatx, vg, vq, theta, cmptcovm[i, ], v = cmptv[i])
    ## lapply(1 :m, bmatx, vg, vq, theta, cmptcovm[i,], v= cmptv[i])
    ## invA <- try(ginv(mA))
    ## if(class(invA) == "try-error"){
    ##     browser()
    ## }
    ## a <- t(invA %*% t(mb))
    geta <- function(i){
       R <- rbind(resp(1, n), cmptresp[, 1] - cmptresp[i, 1], cmptresp[, 3] - cmptresp[i, 3], cmptcovm - matrix(rep(cmptcovm[i, ], n), ncol = n))
       kerall <- dnorm((log(cmptresp[, 1]) - log(cmptresp[i, 1])) / ht1) * dnorm((log(cmptresp[, 3]) - log(cmptresp[i, 3])) / ht2) * apply(dnorm((cmptcovm - matrix(rep(cmptcovm[i, ], n), ncol = n))/hx), 1, prod)
       K  <-   diag(kerall)
       ginv(R %*% K %*% t(R)) %*% R %*% K %*% score [1, ]
   }
    score <- score(cmptresp[, c(1, 3)], theta, cmptcovm)
    a<- do.call(rbind, lapply(1: n, geta))
    apply(score - a, 2, mean)
    #pjscore <-  projscore(vg, vq, theta, cmptresp[i,c("y1", "y2")], orgresp[i, c("y1", "y2")],  cmptcovm[i, ], a,  v = cmptv[i])
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
            missscore <- try(apply( diag(as.numeric(kerdis(x, cmptcovm[ix,  -1, drop = F], hx)), nr, nr) %*% diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( kerdis(x, cmptcovm[ix, -1, drop = F], hx) * cendis[ix] ),  1e-200))
            #missscore <- try(apply( diag(as.numeric(kerx(x, cmptcovm[ix,  , drop = F], hx)), nr, nr) %*% diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( kerx(x, cmptcovm[ix, , drop = F], hx) * cendis[ix] ),  1e-200))
            #missscore <- try(apply(diag(as.numeric(kert(y1, cmptresp[ix, "y1"],  ht)), nr, nr)  %*% diag(cendis[ix], nr, nr) %*%cmptscore[ix, , drop  = F ] , 2, sum) / max(sum( kert(y1, cmptresp[ix, "y1"], ht)  * cendis[ix] ),  1e-200))
            
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
            missscore <- try(apply( diag(as.numeric(kerdis(x, cmptcovm[ix, -1, drop = F ], hx)), nr, nr) %*% diag(cendis[ix], nr, nr)%*%  cmptscore[ix, , drop = F] , 2, sum) / (max(sum(   kerdis(x, cmptcovm[ix, -1, drop =  F], hx) * cendis[ix]), 1e-200) ))
            #missscore <- try(apply( diag(cendis[ix], nr, nr)%*%  cmptscore[ix, , drop = F] , 2, sum) / (max(sum(    cendis[ix]), 1e-200) ))
            }
        else
            missscore <- rep(0, q)
    }
    #if(class(missscore) == "try-error"){
     #   browser()
    #}
    return(missscore)
}
estm1 <- function(theta, resp, survData, covm, n, p, mv = rep(1e-5, n)){
    print(theta)
    theta <- c(abs(theta[1 : 3]) , theta[4:q])
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
    ma[[1]] <<-  creata(theta, cmptresp,  p, mx[1],   cmptv)
    ma[[2]] <<-  creata(theta, cmptresp,  p, mx[2],   cmptv)
    cmptscore <- do.call(rbind, lapply(1 : cn,completescore, theta, cmptresp, survData,  cn, p, ma,  cmptcovm, cmptv))
#    browser()
    if(mn  > 0){
#        browser()
        temp <- (summary(survfit(Surv(survData[, "y2"], 1 - survData[, "d2"] )~1), times= survData[, "y2"], extend=TRUE)$surv)
        cendis <- pmax(temp, min(temp[temp>0]))^(-1)
        missscore <- do.call(rbind, lapply(1 : mn, missingscore, theta, missresp, cmptresp, mn,  p, misscovm, cmptcovm, cmptscore, cendis[cmptix],   missv))
    #browser()
        score <- apply(rbind(cmptscore, missscore), 2, sum)
        }else{
            score <- apply((cmptscore), 2, sum)
        }
    score <- sum((score/n)^2)
    
    
}
estm <- function(theta, resp, survData, covm, n, p, mv = rep(1e-5, n)){
    
    #print(theta)
    theta <- c(abs(theta[1 : 3]) , theta[4:q])
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
    ma[[1]] <<-  creata(theta, cmptresp,  p, mx[1],   cmptv)
    ma[[2]] <<-  creata(theta, cmptresp,  p, mx[2],   cmptv)
    cmptscore <- do.call(rbind, lapply(1 : cn,completescore, theta, cmptresp, survData,  cn, p, ma,  cmptcovm, cmptv))
#    browser()
    if(mn  > 0){
#        browser()
#        temp <- (summary(survfit(Surv(survData[, "y2"], 1 - survData[, "d2"] )~1), times= survData[, "y2"], extend=TRUE)$surv)
        cendis <<- (1 - pexp(survData[, "y2"], 1/cr))^(-1)#pmax(temp, min(temp[temp>0]))^(-1)#temp^(-1)#
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
        r1 <- rbinom(1, 1, 0.6)# x[2] < dg
        nm1 <- rnorm(1, 0, 1)
        p1 <- pnorm(r1 + m1, 0, 1)
    
        covm <-  matrix(c(1,  rbinom(1, 1, p1)), p, 1 )#matrix(1, p, 1)#
    }
    kappa1 <- (abs(theta[1]) )
    kappa2 <- (abs(theta[2]) )
    kappa3 <- (abs(theta[3]) )
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- covm
     
    
    
    g <- r1 * rgamma(1, 1/0.5, 1/0.5) + (1 - r1) * (rgamma(1, 1/6, 1/6) + 2)#(1 - x[2]) * rlnorm(1, 0, 0.75) + x[2] * rweibull(1, 2, scale = 2)#rgamma(1, 1/nu1, scale = nu1)  
   
    lb1 <- g * exp((- t(beta1)%*%x)/kappa1)
    lb2 <- g * exp((- t(beta2)%*%x)/kappa2)
    lb3 <- g * exp((- t(beta3)%*%x)/kappa3)
    a1 <- 1/kappa1
    a2 <- 1/kappa2
    a3 <- 1/kappa3
    p1g2 <- lb2 /(lb1 + lb2)
    r <- rbinom(1,  1, p1g2)
    c <- rexp(1, 1/cr)#runif(1, cen1, cen2)
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
    if(sum(is.na(simdata) )!= 0){
        browser()
    }
    return(c(simdata, covm, g, p1g2))
}

simuRsk2 <- function(i, n, p, nu,  theta,  cen1, cen2 ,covm = NULL){
    
    if(is.null(covm)){
        r1 <- rbinom(1, 1, 0.6)# x[2] < dg
        nm1 <- rnorm(1, 0, 1)
        p1 <- pnorm(r1 + m1, 0, 1)
    
        covm <-  matrix(c(1,  rbinom(1, 1, p1)), p, 1 )#matrix(1, p, 1)#
    }
    kappa1 <- (abs(theta[1]) )
    kappa2 <- (abs(theta[2]) )
    kappa3 <- (abs(theta[3]) )
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- matrix(covm[i, ], ncol = p)
#    g <- rgamma(1, 1/nu1, scale = nu1)  
    lb1 <-  exp((- t(beta1)%*%x)/kappa1) #* 1/kappa1
    lb2 <-  exp((- t(beta2)%*%x)/kappa2) #* 1/kappa2
    lb3 <-  exp((- t(beta3)%*%x)/kappa3) #* 1/kappa3
    a1 <- 1/kappa1
    a2 <- 1/kappa2
    a3 <- 1/kappa3
    p1g2 <- lb2 /(lb1 + lb2)
    r <- rbinom(1,  1, p1g2)
    c <- rexp(1, 1/cr)#runif(1, cen1, cen2)
    if(r == 1){
        u <- runif(1)
        t2 <-   ((u^(-nu ) - 1) / (nu * (lb1 + lb2)))^(1/a2)
        t1 = t2 + 3
        
    }else{
        u <- runif(1)
        t1 <-   ((u^(-nu ) - 1) / (nu * (lb1 + lb2)))^(1/a1)
        u <- runif(1)
        t2 <-  ((u^(- nu / (1 + nu)) * ( 1 + nu * lb1 * t1 ^ a1 + nu * lb2 * t1^a2) - ( 1 + nu * lb1 * t1 ^ a1 + nu * lb2 * t1^a2 - nu * lb3 * t1 ^ a3))/ (nu * lb3))^ (1/a3)
    }
    
    return(c(t1, t2, x))
}

simuRsk3 <- function( n, p, nu,  theta,  cen1, cen2 ,covm = NULL){
    t12 <- matrix(NA, n, 2)
   
    
    
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
        
    
    u <- runif(n - n1)
    t12[!ix, 1] <-   ((u^(-nu ) - 1) / (nu * (lb1[!ix] + lb2[!ix])))^(1/a1)
    u <- runif(n - n1)
    t12[!ix, 2] <-  ((u^(- nu / (1 + nu)) * ( 1 + nu * lb1[!ix] * t12[!ix, 1] ^ a1 + nu * lb2[!ix] * t12[!ix, 1]^a2) - ( 1 + nu * lb1[!ix] * t12[!ix, 1] ^ a1 + nu * lb2[!ix] * t12[!ix, 1]^a2 - nu * lb3[!ix] * t12[!ix, 1] ^ a3))/ (nu * lb3[!ix]))^ (1/a3)
    
    
    return(t12)
}

simuRsk1 <- function(i, n, p, nu,  theta,  cen1, cen2 ,covm = NULL){
    if(is.null(covm)){
        r1 <- rbinom(1, 1, 0.6)# x[2] < dg
        nm1 <- rnorm(1, 0, 1)
        p1 <- pnorm(r1 + m1, 0, 1)
    
        covm <-  matrix(c(1,  rbinom(1, 1, p1)), p, 1 )#matrix(1, p, 1)#
    }
    kappa1 <- (abs(theta[1]) )
    kappa2 <- (abs(theta[2]) )
    kappa3 <- (abs(theta[3]) )
    beta1 <- theta[4 : (3 + p)]
    beta2 <- theta[(4 + p) : (3 + 2 * p)]
    beta3 <- theta[(4 + 2* p) : (3 + 3 * p)]
    x <- covm
#    g <- rgamma(1, 1/nu1, scale = nu1)  
    lb1 <-  exp((- t(beta1)%*%x)/kappa1) #* 1/kappa1
    lb2 <-  exp((- t(beta2)%*%x)/kappa2) #* 1/kappa2
    lb3 <-  exp((- t(beta3)%*%x)/kappa3) #* 1/kappa3
    a1 <- 1/kappa1
    a2 <- 1/kappa2
    a3 <- 1/kappa3
    p1g2 <- lb2 /(lb1 + lb2)
    r <- rbinom(1,  1, p1g2)
    c <- rexp(1, 1/cr)#runif(1, cen1, cen2)
    if(r == 1){
        u <- runif(1)
        t2 <-   ((u^(-nu ) - 1) / (nu * (lb1 + lb2)))^(1/a2)
        t1 = t2 + 3
        
    }else{
        u <- runif(1)
        t1 <-   ((u^(-nu ) - 1) / (nu * (lb1 + lb2)))^(1/a1)
        u <- runif(1)
        t2 <-  ((u^(- nu / (1 + nu)) * ( 1 + nu * lb1 * t1 ^ a1 + nu * lb2 * t1^a2) - ( 1 + nu * lb1 * t1 ^ a1 + nu * lb2 * t1^a2 - nu * lb3 * t1 ^ a3))/ (nu * lb3))^ (1/a3)
    }
    y2 = min(t2, c)
    y1 = min(t1, y2)
    d1 <- as.numeric(t1 < y2)
    d2 <- as.numeric(y2 < c)
    simdata <- cbind(y1, d1, y2, d2)
    colnames(simdata) <- c("y1", "d1", "y2", "d2")
    return(c(simdata, covm,  p1g2))
}

kert <- function(t1, vt2, h){
    dnorm((log(t1) - log(vt2))/h)
}
kerx <- function(x1, vx2, h){
    vx2 <- matrix(vx2, ncol = p-1 )
    x1 <- matrix(rep(x1, nrow(vx2)), ncol = p-1, byrow = T)
    h <- matrix(rep(h, nrow(vx2)), ncol = p-1, byrow = T)
    apply(dnorm((vx2 - x1)/h), 1, prod)
}
kerdis <- function(x1, vx2, h){
    vx2 <- matrix(vx2, ncol = p-1 )
    x1 <- matrix(rep(x1, nrow(vx2)), ncol = p-1, byrow = T)
    h <- matrix(rep(h, nrow(vx2)), ncol = p-1, byrow = T)
    abs(1 - (x1 - vx2))
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
theta <- c(0.5, 0.5, 0.5, -0.5,  -1.1,   -0.2, -1.1,   -0.2, -1.1)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 100

p <- 2
ij <- as.matrix(expand.grid(1 : m, 1 : m))
nu1 <- 0.5
nu <- 0.5
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 1500
up = 20
mx <- matrix(c(0, 1), ncol = p)
dg <- dgamma(0.6, 0.3, 0.3)
cr <- 1000
#survData <- do.call(rbind, lapply(1:n, simuRsk, n, p,  theta, 1, 3))
#survData0 <- do.call(rbind, lapply(1:n, simuRsk1, n, p,nu,  theta0, 300, 400))
                                        #lsurvData <- lapply(1 : 100, simall, 0.8, 1.35)
ma <- vector("list")
sRoot <- function(itr){
    survData <- lsurvData[[itr]]
    
    
    weight <<- rep(1/ng, ng)
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    p1 <- mean(covm[, 2])
    covm1 <<- matrix(cbind(rep(1, ng), rbinom(ng, 1, p1)), ncol = p)#matrix(1, ng, p)#
    XY <<- simuRsk3(ng, p, nu, theta, 300, 400, covm1)#do.call(rbind, lapply(1 :ng, simuRsk2, ng,  p, nu,  theta,  300, 400 , covm1))
    dnom <<- likelihood2(XY, covm1, theta)
    ht <<- sum(resp[, "d1"] == 1)^(-2/15) * bw.nrd0(log(resp[resp[, "d1"] == 1, "y1"]))
    hx <<- n ^ (-2/15) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    vg <<- seq(min(lsurvData[[itr]][, 5 + p]), quantile(lsurvData[[itr]][, 5 + p], 0.95), length.out = m+1)
    

    res <- try(dfsane(theta, estm, method = 3, control = list(tol = 1.e-5, noimp = 20, maxit = 200), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)))
    print(res$convergence)
    return(c(res$par, res$residual))
}
getsimdata <- function(itr){
        
            apply(lsurvData[[itr]][, c(2, 4)], 2, sum)/750
        }
#eventrate <- apply(do.call(rbind, lapply(1:100, getsimdata)), 2, mean)
tsRoot <- function(itr) try(sRoot(itr))
#spg(theta, estm1, gr = NULL,  project = NULL, lower = c(rep(0, 3), rep(-Inf, 3 * p)), upper = rep(Inf, 3 * p), method = 3, projectArgs= NULL, control = list(ftol = 1.e-3 ), quiet = FALSE, resp,survData[, d 1:4],  covm, n, p, rep(min(resp[, 1] /2), n))$par
#res <- mclapply(1:10, tsRoot, mc.cores = 15)

#multiroot(estm, c(rep(1, 6), rep(-0.5, 3)), maxiter = 100,  rtol = 1e-6, atol = 1e-8, ctol = 1e-8,useFortran = TRUE, positive = FALSE,jacfunc = NULL, jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,resp,survData[, 1:4],  covm, n, p)

    
    
delike1 <- function(y1, y2, d1, b1x, b2x, b3x,  g, kappa1, kappa2, kappa3 ){
    dd1x <- -1/kappa1^2 * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-1 / kappa1^2)
    dd2x <- -1/kappa2^2 * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-1 / kappa2^2)
    dd3x <- -1/kappa3^2 * (d1) + ((- exp( (log(y2) - log(g) - b3x) / kappa3^2)) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2))) * (-1/kappa3^2)
    dk1 <- (-2 * (log(y1)-b1x -log(g))/kappa1^3 + (-2/kappa1)) * d1 + (- exp( (log(y1) - log(g) - b1x) / kappa1^2)) * (-2 * (log(y1)-b1x -log(g))/kappa1^3)
    dk2 <- (-2 * (log(y1)-b2x -log(g))/kappa2^3 + (-2/kappa2)) * (1 - d1) + (- exp( (log(y1) - log(g) - b2x) / kappa2^2)) * (-2 * (log(y1)-b2x -log(g))/kappa2^3)
    dk3 <- (-2 * (log(y2)-b3x -log(g))/kappa3^3 + (-2/kappa3)) * d1 + (- exp( (log(y2) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y2)-b3x -log(g))/kappa3^2) - (- exp( (log(y1) - log(g) - b3x) / kappa3^2)) * (-2 * (log(y1)-b3x -log(g))/kappa3^2)
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
    d1 <- as.numeric(vt[, 1] < vt[, 2])
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
#vestm <- jacobian(estm2, x = theta,  method="Richardson", method.args=list(), resp, survData, covm, n, p)
#cr <- 1000
#n <- 500
#lsurvData <- mclapply(1 : 1000, simall,0.3, 1.35, mc.cores = 15)
evalestm <- function(itr){
   #ix <- 1 :n#sample(1:n, n, replace = T)
    survData <- lsurvData[[itr]]#[ix, ]
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    #ht <- n^(-2/15) * bw.nrd0(resp[resp[, "d1"] == 1, "y1"])
    #hx <- n ^ (-2/15) * apply(covm[, , drop = F], 2, bw.nrd0)
#    vg <- quantile(survData[, 5 + p], seq(0, 1, length.out = m))
 #   vq <- dlnorm(vg, 0, 1) /sum(dlnorm(vg, 0, 1) )
    dfsane(theta1, estm2, method = 2, control = list(tol = 1.e-7, noimp = 100 ), quiet = FALSE, resp, survData, covm,  n, p)$par
}

                                        #res <- mclapply(1 : 1000, evalestm, mc.cores = 15)
#res100n <- do.call(rbind, res)
#res1002 <- do.call(rbind, res)
#res250n <- do.call(rbind, res)
#res2502 <- do.call(rbind, res)
#res500n <- do.call(rbind, res)
#res5002 <- do.call(rbind, res)
theta1 <- c(theta, 0.5)
#estm2(theta1, resp, survData[, 1:4], covm, n, p)
## res100 <- res1002[, 1:q]
## res250 <- res2502[, 1:q]
## mres100 <- round(apply(res100, 2, median), 3)
## mres250 <- round(apply(res250, 2, median), 3)
##mres500 <- abs(round(apply(res500, 2, median)[1:q] - theta, 3))
## sdres100 <- round(apply(res100, 2, mad), 3)
## sdres250 <- round(apply(res250, 2, mad), 3)
##sdres500 <- round(apply(res500, 2, mad), 3)
## msres100 <- round((mres100- theta)^2 + sdres100^2, 4)
## msres250 <- round((mres250- theta)^2 + sdres250^2, 4)
#msres500 <- round((mres500)^2 + sdres500[1:q]^2, 4)
## paste(mres100, sdres100, msres100, mres250, sdres250, msres250, sep = "&")
