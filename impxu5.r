
                                        #vl1 is a m x 2 dimensional vector, vl1[, "lambda value"], vl1[, "time"]
library("survival")
library("BB")
fA1 <- function(i, beta, vl, resp, cov){
    ix <- vl[, 2] <= resp[i, "y1"] 
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(beta) %*% cov[i, ])
    else
        A <- 0
    return(A)
}
fA2 <- function(i, beta, vl, resp, cov){
    ix <- vl[, 2] <= resp[i, "y2"] 
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(beta) %*% cov[i, ])
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov){
    ix1 <- vl[, 2] > resp[i, "y1"] & vl[, 2] <= resp[i, "y2"] 
    
    A3 <- 0
    if(sum(ix1) != 0)
        A3 <- sum(vl[ix1, 1])
    
    A3 <- A3  * exp(t(beta) %*% cov[i, ])
}
fB <- function(i, theta, resp){
    1/theta + resp[i, "d1"] + resp[i, "d2"]
    
}




slvl2 <- function( vl1, vl2,  resp,  n, zVec){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    
    fvl <- function(j, vl,  flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 0.0001)
        }else 
            ix <- (y2 >= vl[j, 2]-0.0001)
        
        a <- (sum(zVec[ix]))
        
        if(a != 0 ){
            1/a
        }else{
            browser()
            0
        }
        
        
    }
    
    svl1 <- sapply(1 : m, fvl, vl1,  1)
    svl2 <- sapply(1 : f, fvl, vl2, 2)
                                        # svl3 <- sapply(1: g, fvl, vl3,  0)
    
    return(c(svl1, svl2))
    
}


##################################
##data processing and initializing. Divide data to two set, one containing only the covariates. 
##################################
                                        #resp <- matrix(NA, n, )

################################

####without covariates


mtbb <- 1
nitr <- 10
for(simitr in 1 : nsim){
    print(simitr)
    resp <- lsimresp1[[simitr]]
    B1 <- resp[, 1] * resp[, 2]
    B2 <- resp[, 1] + resp[, 2]   
    colnames(resp) <- c("d1", "d2", "y1", "y2")
    d1 <- resp[, 1]
    d2 <- resp[, 2]
    y1 <- resp[, 3]
    y2 <- resp[, 4]
    p <- 6
    ind1 <- (d1 == 1 )
    
    
    surv1 <- summary(survfit(Surv(resp[, "y1"], resp[, "d1"]) ~ 1))
    surv2 <- summary(survfit(Surv(resp[, "y2"], resp[, "d2"]) ~ 1))
    
    m <- length(surv1$time)
    f <- length(surv2$time)
    
    vl10 <- matrix(NA, m, 2)
    vl20 <- matrix(NA, f, 2)
    
    vl10[, 1]  <- (1 / surv1$n.risk)
    vl20[, 1]  <- (1 / surv2$n.risk)
    
    vl10[, 2] <- surv1$time
    vl20[, 2] <- surv2$time
    
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
                                        # vl3 <- vl30
    theta <- ospg <- spg(0.1, score, gr = NULL, method = 3, project = NULL, lower = 0.01, upper = 5, projectArgs = NULL, control = list(), quiet = FALSE)$par

    print(theta)
    mtbb <- c(mtbb
              , theta)
}
save(theta,mtbb,  file = paste(2, nsim, sep = "_"))
score <- function(t){
    
    for(itr in 1: nitr){
        paras <- sapply(1 : n, getA,  vl1, vl2)
        zVec <- sapply(1:n, postz, t, paras, resp)
        subvl <- slvl2( vl1, vl2,  resp,  n, zVec)###dfsane(subvl, slvl, method = 2, control = list(tol = 1e-5, maxit = 10000, triter = 100), quiet = FALSE, (theta), beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)$par#
        vl1 <- cbind((subvl[1 : m]), vl1[, 2])
        vl2 <- cbind((subvl[(m + 1) : (m + f)]),  vl2[, 2])
        
    }
    paras <- sapply(1 : n, getA,   vl1, vl2)
    
    crit <- -margpartial(t, subvl, paras)
}     

a <- sapply(grid, score)
b <- apply(a^2, 2, sum)

grid <- seq(1.75,  2.1, 0.02)
ospg <- spg(0.1, score, gr = NULL, method = 3, project = NULL, lower = 0.01, upper = 5, projectArgs = NULL, control = list(), quiet = FALSE)



margpartial <- function(theta, vl, paras){
    
    B <- 1/theta + B2
    sum(B1) * log(theta + 1)  - sum(B * log(1 + theta* paras))+ sum(log(vl))#
}




getA <- function(i,   vl1, vl2){
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA1(i, beta1, vl1, resp, cov)
    A2 <- fA2(i, beta2, vl2, resp, cov)
    A <- A1 + A2 
}
postlogz <- function(i, theta, paras, resp){
    digamma(1/theta + resp[i, 1]+ resp[i, 2]) - log(1/theta + paras[i])
}
postz <- function(i, theta, paras, resp){
    (1/theta + resp[i, 1]+ resp[i, 2] )/ (1/theta + paras[i])
}


slvtheta <- function(t, dilogZ){
    dilogZ - n* ( digamma(1/t) ) - n * (log(t) - 1)
}
thetascore <- function(i, theta, paras, resp){
    resp[i, 1] * resp[i, 2] / (theta + 1 ) + 1/(theta^2) * log(1 + theta * paras[i]) - (1 / theta + resp[i, 1] + resp[i, 2]) * paras[i] / ( 1 + theta * paras[i])
}
slvtheU <- function(theta, paras, resp){
    sum( resp[, 1] * resp[, 2] / ( (theta + 1)) + 1/(theta^2) * log(1 + theta * paras) - (1 / theta + resp[, 1] + resp[, 2]) * paras / ( 1 + theta * paras))
}
gammalike <- function(t, resp, paras, zVec){
    sum(-log(dgamma(zVec, 1/t , 1/t )))
}





############simulate samples
simufun <- function(i, t, l1, l2, l3){
    g <- rgamma(1, shape = 1/t, rate = 1/t)
    r <- rbinom(1, 1, 0.5)
    t1 <- rexp(1, g * l1)
    t2 <- rexp(1, g * l2)
    u <- runif(1)
    if(t1 < t2){
        t2 <- -log(exp(-l3 * t1 * g) - (u)* exp(-l3 * t1*g))/(l3*g)
    }else{
        t1 <- -log(exp(-l1 * t2 * g) - (u)* exp(-l1 * t2*g))/(l1*g)
    }
    c <- runif(1, 1, 3)
    c(t1, t2, c)
}
l1 <- 1
l2 <- 1
l3 <- l2

simdata <- t(sapply(1 : n, simufun, 2, l1, l2, l3))
d1 <- simdata[, 1] < simdata[, 2]&simdata[, 1] < simdata[, 3]
d2 <- simdata[, 2] < simdata[, 3]
y2 <- pmin(simdata[, 2], simdata[, 3])
y1 <- pmin(simdata[, 1], y2)
simresp <- cbind(d1, d2, y1, y2)
nsimresp <- simdata
colnames(simresp) <- c("d1", "d2", "y1", "y2")
resp <- simresp

 lsimresp1 <- lsimnresp1 <- vector("list")
for(itr in 1 : nsim){
    simdata1 <- t(sapply(1 : n, simufun, 2, l1, l2, l3))
    d1 <- simdata1[, 1] < simdata1[, 2]&simdata1[, 1] < simdata1[, 3]
    d2 <- simdata1[, 2] < simdata1[, 3]
    y2 <- pmin(simdata1[, 2], simdata1[, 3])
    y1 <- pmin(simdata1[, 1], y2)
    simresp1 <- cbind(d1, d2, y1, y2)
    nsimresp1 <- simdata1
    colnames(simresp1) <- c("d1", "d2", "y1", "y2")

    lsimresp1[[itr]] <- simresp1
    lsimnresp1[[itr]] <- simdata1
    
}







testfU <- function(i, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n,  flag){
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- l1* resp[i, "y1"]#fA12(i, beta1, vl1, resp, cov)
    A2 <- l2 * resp[i, "y1"]#fA12(i, beta2, vl2, resp, cov)
    A3 <- l3 * (resp[i, "y2"] - resp[i, "y1"])#fA3(i, beta3, vl3, resp, cov)

    A <- A1 + A2 + A3
                                        #x <- cov[i, ]
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
testwrapfun <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 2){
    tbb <- c((theta), beta1, beta2, beta3)
    testslvtbb(tbb, vl1, vl2, vl3, resp, cov, n, flag = 2) 
}

testslvl2 <- function( theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    
    
    
    lres <- lapply(1 : n, testfU, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)
    mu <- do.call(rbind, lres)
    fvl <- function(j, vl, submu, flg){
        if(flg == 1){
            ix <- y1 >= vl[j, 2]
        }else 
            ix <- (y2 >= vl[j, 2]) & (y1 < vl[j, 2])
        1 /  (sum(submu[ix]) + 0.0001)
    }
    svl1 <- sapply(1 : m, fvl, vl1, mu[, 1], 1)
    svl2 <- sapply(1 : f, fvl, vl2, mu[, 2], 1)
    svl3 <- sapply(1: g, fvl, vl3, mu[, 3], 0)
    return(c(svl1, svl2, svl3))
    
}

res <- rep(NA, nsim)
for(itr  in 1 : nsim){
    res[itr] <- uniroot(testwrapfun, c(0.0001, 100), beta10, beta20, beta30, vl1, vl2, vl3, lsimresp1[[itr]], cov,n,  2)$root
}













slvl3 <- function(slv, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n,  flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    vl1[, 1] <- exp(slv[1:m])
    vl2[, 1] <- exp(slv[(m+1) : (m + f)])
    vl3[, 1] <- exp(slv[(m + f + 1) : (m + f + g)])
    theta <- exp(slv[m + f +  g+ 1])
    paras <- sapply(1 : n, getA,   beta1, beta2, beta3, vl1, vl2, vl3)
    zVec <- sapply(1:n, postz, theta, paras, resp)
    logZ <- sapply(1:n, postlogz, theta, paras, resp)
    dilogZ <- sum(logZ) - sum(zVec)
    
    fvl <- function(j, vl,  flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 0.0001)
        }else 
            ix <- (y2 >= vl[j, 2] - 0.0001) & (y1 < vl[j, 2])
        
        a <- (sum(zVec[ix]))
        
        if(a != 0 ){
            1/ vl[j, 1] - a
        }else{
            browser()
            0
        }
        
        
    }
    
    
    svl1 <- sapply(1 : m, fvl, vl1,  1) 
    svl2 <- sapply(1 : f, fvl, vl2, 1) 
    svl3 <- sapply(1: g, fvl, vl3,  0) 
    svl4 <- slvtheU(theta, paras, resp)  #slvtheta(theta, dilogZ)#
    return(c(svl1, svl2, svl3, svl4))
    
}


slvl4 <- function(slv, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n,  flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    vl1[, 1] <- (slv[1:m])
    vl2[, 1] <- (slv[(m+1) : (m + f)])
    vl3[, 1] <- (slv[(m + f + 1) : (m + f + g)])
    theta <- (slv[m + f +  g+ 1])
    paras <- sapply(1 : n, getA,   beta1, beta2, beta3, vl1, vl2, vl3)
    zVec <- sapply(1:n, postz, theta, paras, resp)
                                        # logZ <- sapply(1:n, postlogz, theta, paras, resp)
                                        # dilogZ <- sum(logZ) - sum(zVec)
    
    fvl <- function(j, vl,  flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 0.0001)
        }else 
            ix <- (y2 >= vl[j, 2] - 0.0001) & (y1 < vl[j, 2])
        
        a <- (sum(zVec[ix]))
        
        if(a != 0 ){
            1/ a
        }else{
            browser()
            0
        }
        
        
    }
    
    
    svl1 <- sapply(1 : m, fvl, vl1,  1) 
    svl2 <- sapply(1 : f, fvl, vl2, 1) 
    svl3 <- sapply(1: g, fvl, vl3,  0)
    maxlike <- -margpartial(theta, c(svl1, svl2, svl3), paras)
                                        #svl4 <- slvtheU(theta, paras, resp)  #slvtheta(theta, dilogZ)#
    return(c(maxlike))
    
}


theta <- 0.5
mtbb <- 1
nitr <- 100
for(simitr in 1 : nsim){
    print(simitr)
    resp <- lsimresp1[[simitr]]
    B1 <- resp[, 1] * resp[, 2]
    B2 <- resp[, 1] + resp[, 2]   
    colnames(resp) <- c("d1", "d2", "y1", "y2")
    d1 <- resp[, 1]
    d2 <- resp[, 2]
    y1 <- resp[, 3]
    y2 <- resp[, 4]
    p <- 6
    ind1 <- (d1 == 1 )
    
    subgix1 <- (d1 == 0)
    respsub1 <- resp[subgix1, ]
    respsub2 <- resp[ind1, ]

    surv1 <- summary(survfit(Surv(resp[, "y1"], ind1) ~ 1))
    surv2 <- summary(survfit(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~ 1))
    surv3 <- summary(survfit(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~ 1))
    m <- length(surv1$time)
    f <- length(surv2$time)
    g <- length(surv3$time)
    vl10 <- matrix(NA, m, 2)
    vl20 <- matrix(NA, f, 2)
    vl30 <- matrix(NA, g, 2)
    vl10[, 1]  <- (1 / surv1$n.risk)
    vl20[, 1]  <- (1 / surv2$n.risk)
    vl30[, 1]  <- (1 / (surv3$n.risk))
    vl10[, 2] <- surv1$time
    vl20[, 2] <- surv2$time
    vl30[, 2] <- surv3$time
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
    vl3 <- vl30

    subvl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]), (theta)) #c(log(vl10[, 1]), log(vl20[, 1]), log(vl30[, 1]), log(theta))
    subvl <- spg(subvl, slvl4, gr = NULL, method = 2, project = NULL, lower = 0.0001, upper = Inf, projectArgs = NULL, control = list(maxit = 3000), quiet = FALSE, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)$par# dfsane(subvl, slvl3, method = 2, control = list(tol = 1e-5, maxit = 3000, triter = 100), quiet = FALSE,  beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)$par
    vl1 <- cbind(exp(subvl[1 : m]), vl1[, 2])
    vl2 <- cbind(exp(subvl[(m + 1) : (m + f)]),  vl2[, 2])
    vl3 <- cbind(exp(subvl[(m + f + 1) : (m + f + g)]),  vl3[, 2])
    theta <- exp(subvl[m + f + g + 1])
    print(theta)
    mtbb <- c(mtbb, theta)
}

