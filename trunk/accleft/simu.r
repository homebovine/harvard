source("accleft3.r")
set.seed(2008)
m = 25
theta <- c(0.5, 0.5, 0.5, -0.5,  -1.2,     -0.3, -1.0,     -0.5, -1.1)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 500

p <- 2
ij <- as.matrix(expand.grid(1 : m, 1 : m))
nu1 <- 0.5
nu <- 0.5
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 1500
up = 20
cen1 <- 1
cen2 <- 1.5
cr <- 1000 ##noncov 5 &90%60%&2.5&80%55% withcov 10 & 90% 4 & 80% 
mx <- matrix(c(0, 1), ncol = p)
##survData <- do.call(rbind, lapply(1:n, simuRsk, n, p,  theta, 1, 3))
##survData0 <- do.call(rbind, lapply(1:n, simuRsk1, n, p,nu,  theta0, 300, 400))

simall <- function(itr, cen1, cen2){
    set.seed(itr + 2014)
    survData <- do.call(rbind, lapply(1:n, simuRsk, n, p,     theta, cen1, cen2))
    filename <- paste("./simdata/simu",  itr, n, p, m, ng,   cr, sep = "_")
    save(survData, file = filename)
    return(survData)
}
lsurvData <- lapply(1 : 1000, simall, cen1, cen2)
