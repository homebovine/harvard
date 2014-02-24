source("sourcefun.r")
mres <- mvl <- vector("list")
###initial values
beta1 <- 1
beta2 <- 1
beta3 <- 0.5
theta <- 0.5
rtime <- 2 #how many iteratives when performing estimation,
#####initial values
estsim <- function(simitr, beta1, beta2, beta3, theta){
    filenames = paste("./simdata/sim", paste(l1, l2, l3, a, b1, b2, b3, sep = ""),  simitr, sep = "_")
    load(filenames)
    print(simitr)
    covmy <<- matrix(covm, n, p)
    resp <<- simresp1
    B1 <<- resp[, 1] * resp[, 2]
    B2 <<- resp[, 1] + resp[, 2]   
    colnames(resp) <<- c("d1", "d2", "y1", "y2")
    d1 <- resp[, 1]
    d2 <- resp[, 2]
    y1 <- resp[, 3]
    y2 <- resp[, 4]
    ind1 <- which(d1 == 1 )
    ind2 <- which((d1 == 0) & (d2 == 1))
    ind3 <- which((d1 ==1) & (d2 == 1))
    p <- ncol(covmy)
    subgix1 <- (d1 == 0)
    respsub1 <- resp[subgix1, ]
    respsub2 <- resp[ind1, ]

    surv1 <- (coxph(Surv(resp[, "y1"], resp[, "d1"]) ~ (covmy)))
    surv2 <- (coxph(Surv(respsub1[, "y2"], respsub1[, "d2"]) ~ covmy[subgix1]))
    surv3 <- (coxph(Surv(respsub2[, "y2"], respsub2[, "d2"]) ~ covmy[ind1, ]))
    o1 <- order(resp[ind1, 3])
    o2 <- order(resp[ind2, 4])
    o3 <- order(resp[ind3, 4])
    cov1 <<- covmy[ind1, , drop = F][o1, , drop = F]
    cov2 <<- covmy[ind2, , drop = F][o2, , drop = F]
    cov3 <<- covmy[ind3, , drop = F][o3, , drop = F]
    m <<- length(ind1)
    f <<- length(ind2)
    g <<- length(ind3)
    vl10 <- matrix(NA, m, 2)
    vl20 <- matrix(NA, f, 2)
    vl30 <- matrix(NA, g, 2)
    bz10 <- basehaz(surv1)
    bz20 <- basehaz(surv2)
    bz30 <- basehaz(surv3)
    ixbz1 <- which(c(bz10[1, 1], diff(bz10[, 1])) != 0)
    ixbz2 <- which(c(bz20[1, 1], diff(bz20[, 1])) != 0)
    ixbz3 <- which(c(bz30[1, 1], diff(bz30[, 1])) != 0)
    bz1 <- bz10[ixbz1, ]
    bz2 <- bz20[ixbz2, ]
    bz3 <- bz30[ixbz3, ]
    bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
    bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
    bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))

    vl10 <- bz1
    vl20 <- bz2
    vl30 <- bz3
    bz1[, 1] <- c(bz1[1, 1], diff(bz1[, 1]))
    bz2[, 1] <- c(bz2[1, 1], diff(bz2[, 1]))
    bz3[, 1] <- c(bz3[1, 1], diff(bz3[, 1]))
    bb <- c(beta1, beta2, beta3, theta)
    n <- nrow(resp)
    vl1 <- vl10
    vl2 <- vl20
    vl3 <- vl30
    vl <- c((vl10[, 1]), (vl20[, 1]), (vl30[, 1]))

    vl1[, 1] <- vl[1 : m]
    vl2[, 1] <- vl[ (m + 1): (m + f)]
    vl3[, 1] <- vl[ (m + f + 1): (m + f+ g)]
    parasall <- sapply(1 : n, getA, theta,   beta1, beta2, beta3, vl1, vl2, vl3, resp, covmy)
    paras <- parasall[4, ]
    crit <- 0
    broot0 <- broot <- c(0, 0, 0, -0.5)
    for(i in 1 : rtime){
        evl <- dfsane(c(log(vl10[, 1]), log(vl20[, 1]), log(vl30[, 1])), scoremaxnobb, method = 2, control = list(trace = FALSE), quiet = FALSE, resp, covmy, beta1, beta2, beta3, theta, vl1, vl2, vl3)$par
        vl <- exp(evl)
        vl1 <- cbind(vl[1 : m], vl1[, 2])
        vl2 <- cbind(vl[ (m + 1): (m + f)], vl2[, 2])
        vl3 <- cbind(vl[ (m + f + 1): (m + f+ g)], vl3[, 2])
        brootn <- multiroot(wrapscoreindv, broot, maxiter = 100, rtol = 1e-06, atol = 1e-08, ctol = 1e-08, useFortran = TRUE, positive = FALSE, jacfunc = NULL,  jactype = "fullint", verbose = FALSE, bandup = 1, banddown = 1,  resp, covmy, n, vl1, vl2, vl3)$root



        broot <- brootn
        beta1 <- brootn[1]
        beta2 <- brootn[2]
        beta3 <- brootn[3]
        theta <- exp(brootn[4])





    }



    print(c(beta1, beta2, beta3, theta))
    mres <- c(beta1, beta2, beta3, theta)
    mvl <- list(vl1, vl2, vl3)

    return(list(mres, mvl))

}
l1 <- 1###weibull shape parametr for lambda1
l2 <- 1 ######weibull shape for 2
l3 <- 2#######weibull shape for 3
a <- 2
b1 <- c(1)
b2 <- c(1)
b3 <- c(0.5)
nsim <- 100 #######simulation
n <- 250 ######sample size
p <- 1 #########number of covariates
theta <- 0.5
estres <- mclapply(1:nsim, estsim, beta1, beta2, beta3, theta,  mc.cores = 8)
res <- t(sapply(1:nsim, getfromlist, estres, 1))
