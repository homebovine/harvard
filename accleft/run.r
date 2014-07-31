source("accleft4.r")
set.seed(2014)
m1 = 25
m <- 25 

theta <- c(0.5, 0.5, 0.5, -0.6,      -0.3,     -0.5)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 100
sta <- 1
end <- sta + 99
p <- 1
ij <- as.matrix(expand.grid(1 : m, 1 : m))
nu1 <- 0.5
nu <- 0.87
cen1 <- 0.5
cen2 <- 1.5
cr <- 5
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 1500
up = 20
mx <- matrix(c(0, 1), ncol = p)
covm1 <<- matrix(1, ng, p)
sRoot <- function(itr){
    filename <- paste("./simdata/simu", itr, n, p, m, ng, cr, sep = "_")
    load(filename)
    #v <- var(survData[, 6])
    #u <- mean(survData[, 6])
    #shape <- u^2 /v
    #scale <- v/u
    
    weight <<- rep(1/ng, ng)
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    
    ht <<- sum(resp[, "d1"] == 1)^(-2/15) * bw.nrd0(log(resp[resp[, "d1"] == 1, "y1"]))
    hx <<- n ^ (-2/15) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    XY <<- simuRsk3(ng, p, nu, theta, cen1, cen2, covm1)
    
    dnom <<- likelihood2(XY, covm1, theta)
    vg <<- seq(quantile(survData[, 5 + p], 0), quantile(survData[, 5 + p], 0.95), length.out = m+1)# quantile(lsurvData[[itr]][, 5 + p], seq(0.001, 0.999, length.out= m + 1))##qlnorm(seq(0.001, 0.999, length.out = m + 1), 0, 1.5)#
    #vq <<- dgamma(vg, shape = shape, scale = scale)
    #dnom <<- apply(sapply(1:m, dm, vg, XY, covm1, theta), 1, sum) + 1e-6
    #numerator  <- lapply(1:m, num, vg, XY, covm1, theta)
    #sumscore <<- Reduce("+", numerator)
    res <- try(dfsane(theta, estm, method = 3, control = list(tol = 1.e-5, noimp = 20, maxit = 300), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)))
    res <- c(res$par, res$residual)
    #save(res, file =paste("./result/res", itr,  n, p,  m, ng,  cr, sep = "_")) 
    #print(res$convergence)
    return(res)
}
tsRoot <- function(itr) try(sRoot(itr))
res <- mclapply(sta:end, tsRoot, mc.cores = 32)
save(res, file = paste("./result/res", sta, n, p, m, ng, cr, sep = "_"))
