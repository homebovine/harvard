source("accleft6.r")
set.seed(2014)
m = 25
m1 = 100
load("initial1")
#load("inibox")
theta <- theta[1 :12]#c(4.30, 1.32, 3.2, 1.09,      -1.92,     1.15)
q <- length(theta) 
mA <- matrix(NA, m, m)
mb <- matrix(NA, q, m)
#n <- 500
sta <- 1
end <- sta + 99
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
cpvar <- 0
male <- 1
## load("realdata")

## survsd <- sd(survData[, c(1, 3)])
## survData <- cbind(survData[, 1:4], rep(1, 621), survData[, 5])
## survData[, c(1, 3)] <- survData[, c(1, 3)]/survsd

getboot <- function(itr){
     ix <- sample(1 : nrow(survData), nrow(survData), replace = T)
     survData <- survData[ix, ]
     #survData[, c(1, 3)] <- survData[, c(1, 3)]# + runif(nrow(survData), 0, 0.001)
     save(survData,   file = paste("simdata/realdata", "binom", itr, sep = "_"))
 }
## lapply(1 : 1000, getboot)
filename <- paste("simdata/realdata", "binom", 1, sep = "_")
load(filename)
n <- nrow(survData)
rm(survData)
sRoot <- function(itr){
    filename <- paste("simdata/realdata", "binom", itr, sep = "_")#paste("./simdata/simu", itr, n, p, m, ng, cr, sep = "_")
    load(filename)
    #n <<- nrow(survData)
    #v <- var(survData[, 6])
    #u <- mean(survData[, 6])
    #shape <- u^2 /v
    #scale <- v/u
    weight <<- rep(1/ng, ng)
    resp <- survData[, 1:4]
    colnames(resp) <- c("y1", "d1", "y2", "d2")
    covm <- matrix(survData[, 5 : (4+ p)], n, p)
    covm1 <<- matrix(cbind(rep(1, ng), rbinom(ng, 1, mean(survData[, 6]))), ncol = p)
    
    XY <<- simuRsk3(ng, p, nu, theta, 300, 400, covm1)#do.call(rbind, lapply(1 :ng, simuRsk2, ng,  p, nu,  theta,  300, 400 , covm1))
    dnom <<- likelihood2(XY, covm1, theta)
    ht <<- sum(resp[, "d1"] == 1)^(-2/15) * bw.nrd0(log(resp[resp[, "d1"] == 1, "y1"]))
    hx <<- n ^ (-2/15) * apply(covm[, -1, drop = F], 2, bw.nrd0)
    vg <<- seq(min(rg), quantile(rg, 1), length.out = m+1)
    

    
    if(cpvar ==0){
        res <- try(dfsane(theta, estm, method = 3, control = list(tol = 1.e-5, noimp = 20, maxit = 200), quiet = FALSE, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)))
    return(c(res$par, res$residual))
    }else if(cpvar == 1){
        A <- ginv(jacobian(estm, theta, method="Richardson", method.args=list(), resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)))
        V <- varestm(theta, resp,survData[,  1:4],  covm, n, p, rep(min(resp[, 1] /2), n)) /n^2
        estv <- A %*% V %*% t(A)
        return(as.vector(estv))
       
    }
    
 
}
tsRoot <- function(itr) try(sRoot(itr))
res <- mclapply(sta:end, tsRoot, mc.cores = 32)
save(res, file = paste("./result/res3", sta, cpvar, n, p, m, ng, "binombox", sep = "_"))
