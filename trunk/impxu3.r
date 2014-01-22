                                        #vl1 is a m x 2 dimensional vector, vl1[, "lambda value"], vl1[, "time"]
library("survival")
library("BB")
fA12 <- function(i, beta, vl, resp, cov){
    ix <- vl[, 2] <= resp[i, "y1"] + 0.0001
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(beta) %*% cov[i, ])
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov){
    ix1 <- vl[, 2] > resp[i, "y1"] & vl[, 2] <= resp[i, "y2"] + 0.0001
    
    A3 <- 0
    if(sum(ix1) != 0)
        A3 <- sum(vl[ix1, 1])
    
    A3 <- A3  * exp(t(beta) %*% cov[i, ])
}
fB <- function(i, theta, resp){
    1/theta + resp[i, "d1"] + resp[i, "d2"]
    
}




slvl2 <- function( theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, zVec, flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
                                        #lres <- zVec# lapply(1 : n, fU, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, zVec, flag = 0)
    mu <- zVec#do.call(rbind, lres)
    fvl <- function(j, vl,  flg){
        if(flg == 1){
            ix <- y1 >= (vl[j, 2] - 0.0001)
        }else 
            ix <- (y2 >= vl[j, 2] - 0.0001) & (y1 < vl[j, 2])
                                        #b <- sum((y2 < vl[j, 2]) & (y1 < vl[j, 2]))
        a <- (sum(zVec[ix]))
        
        if(a != 0 ){
            1/a
        }else{
            browser()
            0
        }
        
        
    }
    
    svl1 <- sapply(1 : m, fvl, vl1,  1)
    svl2 <- sapply(1 : f, fvl, vl2, 1)
    svl3 <- sapply(1: g, fvl, vl3,  0)
    return(c(svl1, svl2, svl3))
    
}


##################################
##data processing and initializing. Divide data to two set, one containing only the covariates. 
##################################
                                        #resp <- matrix(NA, n, )

################################

####without covariates

theta <- 1
mtbb <- 1
nitr <- 100
for(simitr in 1 : nsim){
    print(simitr)
    set.seed <- 2014
    resp <- lsimresp1[[simitr]]
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
                                        #subtime1 <- quantile(vl10[, 2], seq(0.01, 0.9, 0.005), type = 3)
                                        #subtime2 <- quantile(vl20[, 2], seq(0.01, 0.9, 0.005), type = 3)
                                        #subtime3 <- quantile(vl30[, 2], seq(0.01, 0.9, 0.005), type = 3)
                                        #ix1 <- which(vl10[, 2]%in%subtime1)
                                        #ix2 <- which(vl20[, 2]%in%subtime2)
                                        #ix3 <- which(vl30[, 2]%in%subtime3)
                                        #vl10 <- vl10[ix1, ]
                                        #vl20 <- vl20[ix2, ]
                                        #vl30 <- vl30[ix3, ]
                                        #m <- nrow(vl10)
                                        #f <- nrow(vl20)
                                        #g <- nrow(vl30)
    beta1 <-  beta2 <- beta3<- rep(0, p)
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
    vl3 <- vl30

    
                                        score <- function(theta){
    for(itr in 1: nitr){
                                        #  print(itr)
                                        # print(tbb)
         paras <- sapply(1 : n, getA,   beta1, beta2, beta3, vl1, vl2, vl3)
        zVec <- sapply(1:n, postz, theta, paras, resp)
        logZ <- sapply(1:n, postlogz, theta, paras, resp)
                                        # Z1 <- cbind((Z[1 : m]), vl1[, 2])
                                        # Z2 <- cbind((Z[(m + 1) : (m + f)]),  vl2[, 2])
                                        # Z3 <- cbind((Z[(m + f + 1) : (m + f + g)]),  vl3[, 2])
        subvl <- slvl2(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, zVec, flag = 0)#dfsane(subvl, slvl, method = 2, control = list(tol = 1e-5, maxit = 10000, triter = 100), quiet = FALSE, (theta), beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, n, flag = 0)$par#

        vl1 <- cbind((subvl[1 : m]), vl1[, 2])
        vl2 <- cbind((subvl[(m + 1) : (m + f)]),  vl2[, 2])
        vl3 <- cbind((subvl[(m + f + 1) : (m + f + g)]),  vl3[, 2])
       # paras <- sapply(1 : n, getA,   beta1, beta2, beta3, vl1, vl2, vl3)
        zVec <- sapply(1:n, postz, theta, paras, resp)
        logZ <- sapply(1:n, postlogz, theta, paras, resp)
        dilogZ <- sum(logZ) - sum(zVec)
      #  theta <-  uniroot(slvtheta, c(0.01, 10), dilogZ)$root#uniroot(slvtheU, c(0.01, 10), paras, resp)$root##optim(theta, gammalike, gr= NULL, resp, paras, zVec, method = "L-BFGS-B", lower = 0.01, upper = Inf)$par
        beta1 <- beta10#tbb[2 : (p + 1)]
        beta2 <- beta20#tbb[(p + 2) : (2 * p + 1)]
        beta3 <- beta30#tbb[(2 * p + 2) : (3 * p + 1)]
     #  print(theta)

    }
  -log( prod(zVec^(1/theta + resp[, 1] + resp[, 2]) * exp(- zVec * 1/ theta) * 1/theta^(1/theta)/gamma(1/theta)))
    #dilogZ - n* ( digamma(1/theta) ) - n * (log(theta) - 1)
    }
    print(theta)
    mtbb <- c(mtbb, theta)
}

                                        #wrapfun(theta, beta10, beta20, beta30, vl1, vl2, vl3, resp, cov, n, 2)
                                        #   }
theta <- uniroot(score, c(0.001, 100))$root

mtbb <- c(mtbb, tbb)
theta <- tbb

getA <- function(i,  beta1, beta2, beta3, vl1, vl2, vl3){
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov)
    A2 <- fA12(i, beta2, vl2, resp, cov)
    A3 <- fA3(i, beta3, vl3, resp, cov)
    A <- A1 + A2 + A3
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
    resp[i, 1] * resp[i, 2] / (theta) + 1/(theta^2) * log(1 + theta * paras[i]) - (1 / theta + resp[, 1] + resp[, 2]) * paras[i] / ( 1 + theta * paras[i])
    }
slvtheU <- function(theta, paras, resp){
    sum(sapply(1 : n, thetascore, theta, paras, resp))
    }
gammalike <- function(t, resp, paras, zVec){
    sum(-log(dgamma(zVec, 1/t , 1/t )))
}

                                        #Fine's method
nresp <- realresp
ix <- which(is.na(nresp[, 1])&is.na(nresp[, 2]))
ix1 <- which(is.na(nresp[, 1])&!is.na(nresp[, 2]))
ix3 <- which(is.na(nresp[, 2]))
ix4 <- which(is.na(nresp[, 3]))
nresp[ix, 1:2] <- nresp[ix, 3] + 2
nresp[ix1, 1] <- nresp[ix, 2] + 2
nresp[ix3, 2] <- pmax(nresp[ix3, 1], nresp[ix3, 3], na.rm = T) + 2
nresp[ix4, 3] <- pmax(nresp[ix4, 1], nresp[ix4, 2], na.rm = T) + 2
nr <- nrow(nresp)


fineinner <- function(i, resp, n, num){
    
    y1 <- resp[i, 1]
    y2 <- resp[i, 2]
    c <- resp[i, 3]
    
    inner <- function(j,  n, resp){
        y11 <- resp[j, 1]
        y21 <- resp[j, 2]
        c1 <- resp[i, 3]
        ny1 <- min(y1, y11)
        ny2 <- min(y2, y21)
        nc <- min(c, c1)
        D <- (ny1 < ny2) & (ny2 < nc)
        if(num == 1){
            if(!is.na(D)&D == 1){
                S <- min(ny1, ny2, nc)
                R <- min(ny2, nc)
                W <- n/  sum(resp[, 1] >= S & resp[, 2] >= R)
                W * D * (((y1 - y11)*(y2 - y21) >0 ))
            }else{
                0
            }
        }else{
            if(!is.na(D)&D == 1){
                S <- min(ny1, ny2, nc)
                R <- min(ny2, nc)
                W <- n/  sum(resp[, 1] >= S & resp[, 2] >= R)
                W * D * (1 - ((y1 - y11)*(y2 - y21) >0 ))
            }else{
                0 
            }
        }
                                        #- (1 + t)/(2 + t))
        
    }
    sum(sapply((i + 1) : n, inner,  n, resp))
}
fine <- function( resp, n, num){
    sum(sapply(1: (n-1), fineinner,  resp, n, num))
}
fest <- 2
for(itr in 1:nsim){
    nsimresp <- lsimnresp1[[itr]]
    fest <- c(fest, fine( nsimresp, nr, 1) / fine(nsimresp, nr, 0) - 1)
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
    
    ## u <- runif(1)
    ## if(r == 1){
    ##     t1 <- rexp(1, g * l1)
    ##     t2 <- -log(exp(-l3 * t1 * g) - (u)* exp(-l3 * t1*g))/(l3*g)
    ## }else{
    ##     t2 <- rexp(1, g* l2)
    ##     t1 <- -log(exp(-l1 * t2 * g) - (u)* exp(-l1 * t2*g))/(l1*g)
    ## }
    c <- runif(1, 1, 3)
    c(t1, t2, c)
}
l1 <- 1
l2 <- 3
l3 <- 2

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
        1 /  (sum(submu[ix]) + 0.0001)
    }
    svl1 <- sapply(1 : m, fvl, vl1, mu[, 1], 1)
    svl2 <- sapply(1 : f, fvl, vl2, mu[, 2], 1)
    svl3 <- sapply(1: g, fvl, vl3, mu[, 3], 0)
    return(c(svl1, svl2, svl3))
    
}
uniroot(testwrapfun, c(0.0001, 100), beta10, beta20, beta30, vl1, vl2, vl3, simresp, cov,n,  2)$root
res <- rep(NA, nsim)
for(itr  in 1 : nsim){
    res[itr] <- uniroot(testwrapfun, c(0.0001, 100), beta10, beta20, beta30, vl1, vl2, vl3, lsimresp1[[itr]], cov,n,  2)$root
}
