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


l1 <- 1###weibull shape parametr for lambda1
l2 <- 1 ######weibull shape for 2
l3 <- 1#######weibull shape for 3
a <- 2
b1 <- c(1, 1, 1)
b2 <- c(1, 1, 1)
b3 <- c(0.5, 0.5, 0.5)
nsim <- 100 #######simulation
n <- 250 ######sample size
p <- length(b1) #########number of covariates
theta <- 1
lsimresp1 <- lsimnresp1 <- lcovm <-  vector("list")
set.seed(2013)
for(itr in 1 : nsim){
    covm <- matrix(rnorm(p * n, 0, 1), n, p) #covariance matrix 
    simdata1 <- t(sapply(1 : n, simwei2, theta,  l1, l2, l3, b1, b2, b3, a, covm, 2.5, 10))
    d1 <- simdata1[, 1] < simdata1[, 2]&simdata1[, 1] < simdata1[, 3]
    d2 <- simdata1[, 2] < simdata1[, 3]
    y2 <- pmin(simdata1[, 2], simdata1[, 3])
    y1 <- pmin(simdata1[, 1], y2)
    simresp1 <- cbind(simdata1[, 1], d1, y2, d2)
    colnames(simresp1) <- c("t1", "d1", "t2", "d2")

    lsimresp1[[itr]] <- simresp1
   
    lcovm[[itr]] <- cbind(covm)
#    save(covm, simresp1, file = paste("./simdata/sim", paste(l1, l2, l3, a,  sep = ""),  itr, sep = "_"))
}
