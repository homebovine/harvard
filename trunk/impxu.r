#vl1 is a m x 2 dimensional vector, vl1[, "lambda value"], vl1[, "time"]
library("survival")
library("BB")
fA12 <- function(i, beta, vl, resp, cov){
    ix <- vl[, 2] <= resp[i, "y1"]
    if(sum(ix) != 0)
        A <- sum(vl[ix, 1])* exp(t(beta) %*% cov[i, ])
    else
        A <- 0
    return(A)
}
fA3 <- function(i, beta, vl, resp, cov){
    ix1 <- vl[, 2] <= resp[i, "y1"]
    ix2 <- vl[, 2] <= resp[i, "y2"]
    A1 <- A2 <- 0
    if(sum(ix1) != 0)
        A1 <- sum(vl[ix1, 1])
    if(sum(ix2) != 0)
        A2 <- sum(vl[ix2, 1])
    
    A3 <- (A2 - A1) * exp(t(beta) %*% cov[i, ])
}
fB <- function(i, theta, resp){
    1/theta + resp[i, "d1"] + resp[i, "d2"]
    
}
fU <- function(i, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, flag){
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
   # print(c(commd, theta))
    if(flag == 1){
        u1 <- d1 * d2 /(1 + theta) + 1 / theta ^2 * log(1 + theta * A ) - B * A /(1  + theta * A)
        u2 <- d1 * x - B * x * theta * A1 / commd ##a p \times 1 vector
        u3 <- (1 - d1) * d2 * x - B * x * theta * A2 / commd ##a p \times 1 vector
        u4 <- d1 * d2 * x - B * x * theta * A3 / commd ##a p \times 1 vector
        mu <- matrix(c(u1, u2, u3, u4), nrow = 1) #(3 * p + 1)vector
        return(mu)
    }else{
        u5 <- B * theta * exp(t(beta1) %*% x) / commd
        u6 <- B * theta * exp(t(beta2) %*% x) / commd
        u7 <- B * theta * exp(t(beta3) %*% x) / commd
        return(cbind(u5, u6, u7))
    }
    
}
slvtbb <- function(tbb, vl1, vl2, vl3, resp, cov, flag = 1){
    theta <- exp(tbb[1]) 
    beta1 <- tbb[2 : (p + 1)]
    beta2 <- tbb[(p + 2) : (2 * p + 1)]
    beta3 <- tbb[(2 * p + 2) : (3 * p + 1)]
    lres <- lapply(1 : n, fU,theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, flag = 1)
    mu <- do.call(rbind, lres)
    scu <- apply(mu, 2, sum)
}
slvl <- function(subvl, theta, beta1, beta2, beta3, vl1, vl2, vl3, subvl1, subvl2, subvl3, resp, cov, flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    subvl1[, 1] <- exp(subvl[1 : m1])
    subvl2[, 1] <- exp(subvl[(m1 + 1) : (m1 + f1)])
    subvl3[, 1] <- exp(subvl[(m1 + f1 + 1) : (m1 + f1 + g1)])
    vl1[, 1] <- approxfun(subvl1[, 2], subvl1[, 1])(vl1[, 2])
    vl2[, 1] <- approxfun(subvl2[, 2], subvl2[, 1])(vl2[, 2])
    vl3[, 1] <- approxfun(subvl3[, 2], subvl3[, 1])(vl3[, 2])
    lres <- lapply(1 : n, fU, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, flag = 0)
    mu <- do.call(rbind, lres)
    fvl <- function(j, vl, submu, flg){
        if(flg == 1){
            ix <- y1 >= vl[j, 2]
        }else
            ix <- (y2 >= vl[j, 2]) & (y1 < vl[j, 2])
        1 / vl[j, 1] - sum(submu[ix])
    }
    svl1 <- sapply(1 : m1, fvl, subvl1, mu[, 1], 1)
    svl2 <- sapply(1 : f1, fvl, subvl2, mu[, 2], 1)
    svl3 <- sapply(1: g1, fvl, subvl3, mu[, 3], 0)
    return(c(svl1, svl2, svl3))
    
}

slvl1 <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, rvl1, rvl2, rvl3, resp, cov, flag = 0){
    d1 <- resp[, "d1"]
    d2 <- resp[, "d2"]
    y1 <- resp[, "y1"]
    y2 <- resp[, "y2"]
    lres <- lapply(1 : n, fU, theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, flag = 0)
    mu <- do.call(rbind, lres)
    lres <- lapply(1 : n, fU, theta, beta1, beta2, beta3, rvl1, rvl2, rvl3, resp, cov, flag = 0)
    mu1 <- do.call(rbind, lres)
    fvl <- function(j, vl, submu, flg){
        if(flg == 1){
            ix <- y1 >= vl[j, 2]
        }else
            ix <- (y2 >= vl[j, 2]) & (y1 < vl[j, 2])
        1 / vl[j, 1] - sum(submu[ix])
    }
    svl1 <- max(abs(sapply(1 : m, fvl, vl1, mu[, 1], 1) - sapply(1 : m, fvl, rvl1, mu1[, 1], 1)))
    svl2 <- max(abs(sapply(1 : f, fvl, vl2, mu[, 2], 1) - sapply(1 : f, fvl, rvl2, mu1[, 2], 1)))
    svl3 <- max(abs(sapply(1: g, fvl, vl3, mu[, 3], 0) - sapply(1: g, fvl, rvl3, mu1[, 3], 0)))
    return(c(svl1, svl2, svl3))
    
}

##################################
##data processing and initializing. Divide data to two set, one containing only the covariates. 
##################################
#resp <- matrix(NA, n, )

################################
set.seed <- 2014
myData <- read.csv("DataForSebastien121127.csv")
myData <- myData[complete.cases(myData[, 1:9]), ]
n <- nrow(myData)
ix <- sample(1:n, 300)
myData <- myData[ix, ]
y1 <- myData$DiagnosisAge - myData$EnrollAge
y2 <- myData$AgeAtDeath - myData$EnrollAge


d1 <- !is.na(y1) 
d2 <- !is.na(y2)
cen1 <- runif(sum((1 - d2) * d1), min(myData$AgeAtDeath - myData$DiagnosisAge, na.rm = T), max(myData$AgeAtDeath - myData$DiagnosisAge, na.rm = T) + 10)
cen2 <- runif(sum((1 - d2) * (1 - d1)), min(myData$DiagnosisAge - myData$EnrollAge, na.rm = T), max(myData$DiagnosisAge - myData$EnrollAge, na.rm = T) + 10)
#cenAge <- myData$EnrollAge[1 - d2] + cen1
y2[which( d2 == 0&(d1 ==1))] <- cen1 + y1[which( d2 == 0&(d1 ==1))]
y2[which( d2 == 0&(d1 ==0))] <- cen2
y1[which(d1 == 0)] <- y2[which(d1 == 0)]
ixnoz <-y1 != 0&y2 !=0
resp <- cbind(as.numeric(d1), as.numeric(d2), y1, y2)
resp <- resp[ixnoz, ]
cov <- cbind(myData$gender, myData$race, myData$education, myData$marital, myData$CESD_Score, as.numeric(myData$apoe_raw))
cov <- cov[ixnoz, ]
d1<- d1[ixnoz ]
d2 <- d2[ixnoz]
y1<- y1[ixnoz ]
y2 <- y2[ixnoz]

colnames(resp) <- c("d1", "d2", "y1", "y2")
p <- ncol(cov)
ind1 <- (d1 == 1)
m <- sum(ind1)
ind2 <- (d2 == 1 & d1 == 0)
f <- sum(ind2)
ind3 <- (d2 == 1& d1 == 1)
g <- sum(ind3)
vl10 <- matrix(NA, m, 2)
vl20 <- matrix(NA, f, 2)
vl30 <- matrix(NA, g, 2)
m1 <- min(100, m)
f1 <- min(100, f)
g1 <- min(100, g)
subinx1 <- sample(1 : m, m1)
subinx2 <- sample(1 : f, f1)
subinx3 <- sample(1 : g, g1)
vl10[, 2] <- resp[ind1, "y1"]
vl20[, 2] <- resp[ind2, "y2"]
vl30[, 2] <- resp[ind3, "y2"]


surv1 <- summary(survfit(Surv(resp[, "y1"], ind1) ~ 1),  times = resp[ind1, "y1"], extend = TRUE)
surv2 <- summary(survfit(Surv(resp[, "y2"], ind2) ~ 1),  times = resp[ind2, "y2"], extend = TRUE)
surv3 <- summary(survfit(Surv(resp[, "y2"], ind3) ~ 1),  times = resp[ind3, "y2"], extend = TRUE)
vl10[, 1]  <- (1 / surv1$n.risk)
vl20[, 1]  <- (1 / surv2$n.risk)
vl30[, 1]  <- (1 / surv3$n.risk)
subvl10 <- vl10[subinx1, ]
subvl20 <- vl20[subinx2, ]
subvl30 <- vl30[subinx3, ]
theta0 <- 0
beta10 <-  beta20 <- beta30 <- rep(0, p)

n <- nrow(resp)

tbb <- c(theta0, beta10, beta20, beta30)
vl <- c(log(vl10[, 1]), log(vl20[, 1]), log(vl30[, 1]))
subvl <- c(log(subvl10[, 1]), log(subvl20[, 1]), log(subvl30[, 1]))
vl1 <- vl10
vl2 <- vl20
vl3 <- vl30
subvl1 <- cbind((subvl[1 : m1]), resp[subinx1, "y1"])
subvl2 <- cbind((subvl[(m1 + 1) : (m1 + f1)]), resp[subinx2, "y2"])
subvl3 <- cbind((subvl[(m1 + f1 + 1) : (m1 + f1+ g1)]), resp[subinx3, "y2"])
nitr <- 20
mtbb <- tbb
for(itr in 1: nitr){
    print(itr)
    tbb <- dfsane(tbb, slvtbb, method = 2, control = list(tol = 1e-5,  maxit = 1500, triter = 100), quiet = FALSE, vl1, vl2, vl3, resp, cov, flag = 1)$par
    theta <- tbb[1]
    beta1 <- tbb[2 : (p + 1)]
    beta2 <- tbb[(p + 2) : (2 * p + 1)]
    beta3 <- tbb[(2 * p + 2) : (3 * p + 1)]
    subvl <- dfsane(subvl, slvl, method = 2, control = list(tol = 1e-5, maxit = 1500, triter = 100), quiet = FALSE, exp(theta), beta1, beta2, beta3, vl1, vl2, vl3, subvl1, subvl2, subvl3,  resp, cov, flag = 0)$par
    
    subvl1 <- cbind((subvl[1 : m1]), resp[subinx1, "y1"])
    subvl2 <- cbind((subvl[(m1 + 1) : (m1 + f1)]), resp[subinx2, "y2"])
    subvl3 <- cbind((subvl[(m1 + f1 + 1) : (m1 + f1 + g1)]), resp[subinx3, "y2"])
    vl1[, 1] <- approxfun(subvl1[, 2], exp(subvl1[, 1]))(vl1[, 2])
    vl2[, 1] <- approxfun(subvl2[, 2], exp(subvl2[, 1]))(vl2[, 2])
    vl3[, 1] <- approxfun(subvl3[, 2], exp(subvl3[, 1]))(vl3[, 2])
   # vl2 <- approxfunc(subvl2[, 1], subvl2[, 2])(resp[ind2, "y2"])
    #vl1 <- cbind(exp(vl[1 : m]), resp[ind1, "y1"])
   # vl2 <- cbind(exp(vl[(m + 1) : (m + f)]), resp[ind2, "y2"])
    #vl3 <- cbind(exp(vl[(m + f + 1) : (m + f + g)]), resp[ind3, "y2"])
    mtbb <- c(mtbb, tbb)
}

theta <- exp(tbb[1])
beta1 <- tbb[2 : (p + 1)]
beta2 <- tbb[(p + 2) : (2 * p + 1)]
beta3 <- tbb[(2 * p + 2) : (3 * p + 1)]
#vl1 <- cbind(exp(vl[1 : m]), resp[ind1, "y1"])
#vl2 <- cbind(exp(vl[(m + 1) : (m + f)]), resp[ind2, "y2"])
#vl3 <- cbind(exp(vl[(m + f + 1) : (m + f + g)]), resp[ind3, "y2"])


######variance calculation
#####The finite dimensional parameters. matrix I_{11}
ufun <- function(i, resp, cov, theta, beta1, beta2, beta3, vl1, vl2, vl3, flg1, flg2 ){
    d1 <- resp[i, "d1"]
    d2 <- resp[i, "d2"]
    A1 <- fA12(i, beta1, vl1, resp, cov)
    A2 <- fA12(i, beta2, vl2, resp, cov)
    A3 <- fA3(i, beta3, vl3, resp, cov)
     A <- A1 + A2 + A3
    x <- cov[i, ]
    B <- fB(i, theta, resp)
    commd <- (1 + theta * A)
    if(flg1 == 1 & flg2 == 1){
        -d1 * d2 / (1 + theta)^2 - 2/theta^3 *log(1 + theta * A) + 2 / theta^2 * A/(1 + theta * A) + B * (A/(1 + theta * A))^2
    }else if(flg1 == 1 & flg2 == 2){
        (1 /( theta* (1 + theta * A)) - B /(1 + theta * A)^2) *x *A1
    }else if(flg1 == 1 & flg2 == 3){
        (1 /( theta* (1 + theta *  A)) - B /(1 + theta * A)^2) *x *A2
    }else if(flg1 == 1 & flg2 == 4){
        (1 /( theta* (1 + theta * A)) - B /(1 + theta * A)^2) *x *A3
    }else if(flg1 == 2 & flg2 == 2){
       as.numeric( -B * theta * A1 * (1 + theta * A2 + theta * A3)/ (1 + theta * A)^2) * (x%*% t(x))
    }else if(flg1 == 2 & flg2 == 3){
        as.numeric(B * theta * A1 * ( theta * A2)/ (1 + theta * A)^2) *  x%*% t(x)
    }else if(flg1 == 2 & flg2 == 4){
        as.numeric(B * theta * A1 * ( theta * A3)/ (1 + theta * A)^2) *  x%*% t(x)
    }else if(flg1 == 3 & flg2 == 3){
        as.numeric(-B * theta * A2 * ( 1 + theta * A1  + theta * A3)/ (1 + theta * A)^2) * x%*% t(x)
    }else if(flg1 == 3 & flg2 == 4){
        as.numeric(B * theta * A2 * ( theta * A3)/ (1 + theta * A)^2) * x%*% t(x)
    }else if(flg1 == 4 & flg2 == 4){
        as.numeric(-B * theta * A3 * (1 + theta * A1 +  theta * A2)/ (1 + theta * A)^2) * x%*% t(x)
    }
}

ulambda <- function(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, xi){
    rvl1 <- cbind(vl1[, 1] + xi, vl1[, 2])
    rvl2 <- cbind(vl2[, 1] + xi, vl2[, 2])
    rvl3 <- cbind(vl3[, 1] + xi, vl3[, 2])
    tbb <- c(log(theta), beta1, beta2, beta3)
    s1 <- slvtbb(tbb, vl1, vl2, vl3, resp, cov, 1)
    s2 <- slvtbb(tbb, rvl1, vl2, vl3, resp, cov, 1)
    ul1 <- (s2 - s1) / xi
    sl1 <- slvl1(theta, beta1, beta2, beta3, vl1, vl2, vl3,  rvl1, vl2, vl3, resp, cov, flag = 0) / xi
    
    
   
    s2 <- slvtbb(tbb, vl1, rvl2, vl3, resp, cov, 1)
    ul2 <- (s2 - s1) / xi
    sl2 <- slvl1(theta, beta1, beta2, beta3, vl1, vl2, vl3,  vl1, rvl2, vl3,   resp, cov, flag = 0) / xi
    
    
    
    s2 <- slvtbb(tbb, vl1, vl2, rvl3, resp, cov, 1)
    ul3 <- (s2 - s1) / xi
    sl3 <- slvl1(theta, beta1, beta2, beta3, vl1, vl2, vl3,  vl1, vl2, rvl3,   resp, cov, flag = 0) / xi
    msl <- rbind(sl1, sl2, sl3)
    msl[lower.tri(msl)]  <-  msl[upper.tri(msl)]
   # browser()
    finmax <- cbind(rbind(I11, ul1, ul2, ul3), rbind(cbind(ul1, ul2,  ul3), msl))
    return(finmax)
}









xi <- 0.0001

I11 <- matrix(NA, 19, 19)
I11[1, 1] <- sum(sapply(1 : n, ufun, resp, cov, theta, beta1, beta2, beta3, vl1, vl2, vl3, 1, 1))
for (k in 1 : 3){
    I11[1, (p * (k-1) + 2) : (p* (k - 1) + p + 1)] <- I11[(p * (k-1) + 2) : (p * (k - 1) + p + 1), 1] <- apply(do.call(rbind, lapply(1 : n, ufun, resp, cov, theta, beta1, beta2, beta3, vl1, vl2, vl3, 1, k + 1)), 2, sum)
}

for(k in 1 : 3){
    for (l in  (k : 3)){
        I11[((p * (k - 1) + 2): (p* (k - 1) + p + 1)), ((p * (l - 1) + 2): (p* (l - 1) + p + 1))]<- I11[((p * (l - 1) + 2): (p* (l - 1) + p + 1)), ((p * (k - 1) + 2): (p* (k - 1) + p + 1))] <-  Reduce("+", lapply(1 : n, ufun, resp, cov, theta, beta1, beta2, beta3, vl1, vl2, vl3, k + 1, l + 1))
    }
}
Imtx <- ulambda(theta, beta1, beta2, beta3, vl1, vl2, vl3, resp, cov, xi)

coxd1 <- coxph(Surv(resp[, "y1"], resp[, "d1"] )~cov[, 1] + cov[, 2] + cov[, 3] + cov[, 4]+ cov[, 5]+cov[, 6])
coxd2 <- coxph(Surv(resp[, "y2"], (1 - resp[, "d1"]) * resp[, "d2"])~cov[, 1] + cov[, 2] + cov[, 3] + cov[, 4]+ cov[, 5]+cov[, 6])
coxd3 <- coxph(Surv(resp[, "y2"], resp[, "d1"] * resp[, "d2"])~cov[, 1] + cov[, 2] + cov[, 3] + cov[, 4]+ cov[, 5]+cov[, 6])
coxd4 <- coxph(Surv(resp[, "y1"], resp[, "d1"] * (1 - resp[, "d2"]))~cov[, 1] + cov[, 2] + cov[, 3] + cov[, 4]+ cov[, 5]+cov[, 6])
coxd5 <- coxph(Surv(resp[, "y1"], resp[, "d2"] )~cov[, 1] + cov[, 2] + cov[, 3] + cov[, 4]+ cov[, 5]+cov[, 6])
#######reconstruct data
ix <- which(resp[, 1] == 0)
ix1 <- which(resp[, 1] != 0)
l1 <- length(ix)
l2 <- length(ix1)
colnames(resp3) <- c("id", "d", "t", "death")
resp2 <- resp[ix, c(1, 4)]
resp1 <- rbind(resp[ix1, c(1, 3)], cbind(rep(0, l2), resp[ix1, c(4)]))
resp1 <- cbind(c(1 : l2, 1 : l2 ), resp1, c(rep(0, l2), rep(1, l2)))
resp2 <- cbind(c((l2+1) : (l2 + l1)), resp2, rep(1, l1))
resp3 <- rbind(resp1, resp2)
cov1 <- rbind(cov[ix1, ], cov[ix1, ], cov[ix, ])
o <- order(resp3[, "id"])
cov1 <- cov1[o, ]
resp3 <- resp3[o, ]

resp3 <- data.frame(resp3)
a <- frailtyPenal(Surv(t, d) ~ cluster(id) + cov1[, 1] + cov1[, 2] + cov1[, 3] + cov1[, 4]+ cov1[, 5]+cov1[,6] + terminal(death), formula.terminalEvent = ~ cov1[, 1] + cov1[, 2] + cov1[, 3] + cov1[, 4]+ cov1[, 5]+cov1[,6], data = resp3, Frailty = TRUE, joint = TRUE, n.knots=5, kappa1=1000, kappa2=1000)
a <- frailtyPenal(Surv(t, d) ~ + cov1[, 1] + cov1[, 2] + cov1[, 3] +death, data = resp1,  n.knots=5, kappa1=9.55e9, kappa2=1.41e12)





















