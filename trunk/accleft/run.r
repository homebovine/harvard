source("accleft3.r")
set.seed(2014)
m = 25
m1 = 100
theta <- c(0.5, 0.5, 0.5, -0.6,    -0.3,   -0.5)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
n <- 100
sta <- 1
end <- sta + 99
p <- 1
ij <- as.matrix(expand.grid(1 : m, 1 : m))
nu1 <- 0.5
nu <- 0.5
#ij <- ij[ij[, 1] >= ij[, 2], ]
ng <- 1500
up = 20
mx <- matrix(c(0, 1), ncol = p)
sRoot <- function(itr){
    filename <- paste("./simdata/simu", itr, n, p, m, sep = "_")
    load(filename)
    #cx <<- gaussLegendre(ng, quantile(as.vector(survData[, c(1)]), 0), quantile(as.vector(survData[, c(1)]), 1))
    #cx <<- gaussLegendre(ng, 0.01, 20)
    #x <<- cx$x
    #wx <<- cx$w
    #cy <<- gaussLegendre(ng, quantile(as.vector(survData[, c( 3)]), 0), quantile(as.vector(survData[, c( 3)]), 1))
    #cy <<- gaussLegendre(ng, 0.01, 20)
    #y <<- cy$x
    #wy <<- cy$w
    #mgrid <<- meshgrid(x, y)
    #cbind(runif(ng, 0, up), runif(ng, 0, up))
    
    weight <<- rep(1/ng, ng)
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    covm1 <<- matrix(1, ng, p)#matrix(cbind(rep(1, ng), rnorm(ng, 0, 1)), ncol = p)
    XY <<- simuRsk3(ng, p, nu, theta, 300, 400, covm1)#do.call(rbind, lapply(1 :ng, simuRsk2, ng,  p, nu,  theta,  300, 400 , covm1))
    dnom <<- likelihood2(XY, covm1, theta)
    ht <<- sum(resp[, "d1"] == 1)^(-2/15) * bw.nrd0(log(resp[resp[, "d1"] == 1, "y1"]))
    hx <<- n ^ (-2/15) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    vg <<- seq(min(lsurvData[[itr]][, 6]), max(lsurvData[[itr]][, 6]), length.out = m+1)#qgamma(seq(0.000001, 0.99999999999999, length.out = m + 1), 1/nu, 1/nu)
    vg1 <<- qgamma(seq(0.00000001, 0.99999999999, length.out = m1 + 1), 1/nu, 1/nu)
    

    res <- try(dfsane(theta, estm, method = 3, control = list(tol = 1.e-5, noimp = 20, maxit = 200), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n))$par)
    save(res, file = paste("./result/res", itr, n, p, sep = "_"))
    return(res)
}
tsRoot <- function(itr) try(sRoot(itr))
res <- mclapply(sta:end, tsRoot, mc.cores = 15)
save(res, file = paste("./result/res", sta, n, p, m, sep = "_"))
