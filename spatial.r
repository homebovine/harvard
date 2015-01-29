library(dlnm)
inidate <- as.Date("1999-01-01", format = "%Y-%m-%d")
lag7 <- function(i, loc, subdataHW, datahwd, datactrl){
    subdataHW$date <- as.Date(subdataHW$date)
    subdataHW[i, 6:12]<- seq(subdataHW$date[i]-6, subdataHW$date[i], 1)%in% datahwd
    subdataHW[i, 13] <- subdataHW$date[i]%in% datactrl
    subdataHW[i, 6:13]
    
                          

}

uniloc <- unique(dataHWcontrol$FIPS)
len <- length(uniloc)
colnames(dataHW)[5] <- "admit"
loccoef <- function(i){
    load(paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    colnames(subdataHWcontrol)[5] <- "admit"
 
    mQ <- as.matrix((subdataHWcontrol[, seq(12, 6, -1)]), ncol = 7)
    cov1 <- crossbasis(mQ, lag = c(0, 6), argvar = list(fun = "lin", cen = TRUE), arglag = list(fun = "poly", degree = df) )
    cov2 <- log(as.numeric(subdataHWcontrol$date - inidate)/365)
    glmres <- glm(admit~cov1  +cov2 +  offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    glmcoef <- summary(glmres)$coefficients
    C <- solve(t(mQ) %*% mQ) %*% t(mQ) %*% cov1
    
    bb <- sum(abs(C %*% glmcoef[2:(df + 2), 1]))
    vbb <- C %*% vcov(glmres)[2:(df + 2), 2:(df + 2)] %*% t(C)
    sgn <- sign(C %*% glmcoef[2:(df + 2), 1])
    sdbb <-sqrt(t(sgn) %*% vbb %*% (sgn))# sqrt(nrow(subdataHWcontrol))^(-1)#
    rho <- C%*%glmcoef[2:(df + 2), 1]/bb
    
    
    c(bb, sdbb, rho)
}

savedata <- function(i){
    subdataHW <- dataHW[dataHW$FIPS == uniloc[i], ]
    
    datahwd <- dataHWcontrol$date[dataHWcontrol$hw3 &dataHWcontrol$FIPS == uniloc[i]]
    datactrl <- dataHWcontrol$date[dataHWcontrol$ctrl3 &dataHWcontrol$FIPS == uniloc[i]]
    save(subdataHW, datahwd, datactrl, file = paste("./spatial/loc55", uniloc[i], sep = "_"))
    
}

getdata <- function(i){
    load(paste("./spatial/loc55", uniloc[i], sep = "_"))
    colnames(subdataHW)[5] <- "admit"
    subdataHW <- cbind(subdataHW, matrix(NA, nrow(subdataHW), 8))
    subdataHW[, 6:13] <- matrix(do.call(rbind, lapply(1 : nrow(subdataHW), lag7,  uniloc[i], (subdataHW), datahwd, datactrl)), ncol = 8, byrow = T)
    subdataHWcontrol <- subdataHW[apply(subdataHW[, 6:13], 1, sum) > 0, ]

    print(i)
    save(subdataHWcontrol, file = paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    return(subdataHWcontrol)
    
}


trygetdata <- function(i){
    res<- try(getdata(i))
    if(class(res) == "try-error"){
        return(NULL)
        }else{
            return(res)
            }
    
    }
lapply(1:len, savedata)
aggdata <- mclapply(1 : len, trygetdata, mc.cores = 10)
aggallData <- do.call(rbind, aggdata)
meanhwday <- aggregate(aggallData[, 6:12], by=list(aggallData$FIPS), FUN = mean)
sumcount <- aggregate(aggallData[, 5], by=list(aggallData$FIPS), FUN = sum)
sumadjust <- aggregate(aggallData[, 5] / aggallData[, 4], by=list(aggallData$FIPS), FUN = sum)
spsz <- aggregate(aggallData[, 5], by=list(aggallData$FIPS), FUN = length)

tryloc <- function(i){
    res <- try(loccoef(i))
    
    if(class(res) != "try-error"){
        return(res)
    }else{
        return(rep(NA, 2 + 7))
        }
    }
aggbb <- lapply(1 : len, tryloc)#, mc.cores = 10)
aggbb1 <- aggbb <- do.call(rbind, aggbb)

aggpvalue <- mclapply(1 : len, tryloc, mc.cores = 10)
aggpvalue <- do.call(rbind, aggpvalue)



ix <- which(sumcount[, 2] > 60)
uniloc <- uniloc[ix]
len <- length(uniloc)



AltLong <- AltLong[ix, ]
AltLong <- as.matrix(AltLong)
h <- c(bw.nrd0(AltLong[, 1]),bw.nrd0(AltLong[, 2])) * len^(-1/15)
dis <- dis[ix, ix]
newaa <- matrix(1, ncol= 1, nrow= len)
aggbb <- aggbb1[ix, ]
object <- function(aa, s, invSigma, Kh){

 
    estm <- aggbb[, 1] - aa[1]# - (s - AltLong) %*% aa[2:3]#apply((aa[, 2 : 3]) *  (s - AltLong), 1, sum)
    temp <- (t(estm) %*% Kh %*% invSigma %*% estm) /len
    temp
    
}

objectrho <- function(mrho,   a0){
    estm <- aggbb[, 1] - a0
    rho <- mrho[1]
    sigma<- mrho[2]
   
    #mrho <-  1 - 3 * abs(dis)/(2 * rho) + abs(dis)^3/ (2 * rho ^3)
    #mrho[dis > rho] <- 0
     mrho <-     exp(-sqrt(5)  * dis  /  rho) * (1 + sqrt(5) * dis  /rho + 5 * (dis) ^2 /(3 * rho ^2))
        mrho[dis > rho] <- 0
    #mrho <-  (1 + sqrt(3) * dis /rho ) * exp(-sqrt(3) * dis /rho)
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho * sigma
    invSigma <- ginv(Vpara)
   (1/2 *  t(estm)  %*% invSigma %*% estm +1/2 * determinant(Vpara, logarithm = TRUE)[[1]])/len
    
}

objectrho1 <- function(rho,   a0){
    
    
    sigma <- rho[2]
    rho <- rho[1]
    #mrho <- exp(-sqrt(5)  * dis  /  rho) * (1 + sqrt(5) * (dis) /rho + 5 * (dis) ^2 /(3 * rho ^2))
    #mrho <-  (1 + sqrt(3) * dis /rho ) * exp(-sqrt(3) * dis /rho)
    mrho <-  1 - 3 * abs(dis)/(2 * rho) + abs(dis)^3/ (2 * rho^3)
    mrho[dis > rho] <- 0
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho #* sigma
    estm <- (aggbb[, 1]  -  a0)%*%t(aggbb[, 1]  -  a0)   - Vpara
    max(eigen(t(estm)  %*% estm)$values)#+1/2 * log(det(Sigma))
    
}


temp <- optim(aa, object, gr= NULL, AltLong[i, ], method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.1), upper = c(Inf, Inf, Inf, 0.9))
AltLong <- cbind(rep(AltLong[, 1], each = 7), rep(AltLong[, 2], each = 7))
getaa <- function(i, aa, rho, sigma){
 #   print(i)
    s <- AltLong[i, ]
    
    s <- matrix(rep(s, len), ncol = 2, byrow = T)
    Kh <- diag((dnorm((s[, 1] - AltLong[, 1])/ h[1]) * dnorm((s[, 2] - AltLong[, 2])/ h[2])))
    #mrho <-  1 - 3 * abs(dis)/(2 * rho) + abs(dis)^3/ (2 * rho^3)
    #mrho[dis > rho] <- 0

    mrho <-  exp(-sqrt(5)  * dis /  rho) * (1 + sqrt(5) * dis /rho + 5 * (dis) ^2 /(3 * rho ^2))#(1 + sqrt(3) * dis /rho ) * exp(-sqrt(3) * dis /rho)#
        mrho[dis > rho] <- 0
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho * sigma
    invSigma <- ginv(Vpara)
    aa <- aa[i]
    aa <- spg(aa, object, gr = NULL, method = 3, project = NULL, lower = 0.001, upper = Inf, projectArgs = NULL, control = list(trace = TRUE), quiet = FALSE, s, invSigma, Kh )$par
}

estm <- lapply(1:len, getaa, aa, rho1)
tol <- 1e-5
for(j in 1: nsim){
    rho <- newrho[1]
    sigma <- newrho[2]
    
    aa <- newaa
    newaa <- do.call(rbind, mclapply(1:len, getaa, aa, rho, sigma, mc.cores = 20))
    newrho <- spg(newrho, objectrho, gr = NULL, method = 3, project = NULL, lower = c(0.001, 0.001), upper = c(100,  100), projectArgs = NULL, control = list(trace = TRUE), quiet = FALSE, newaa)$par
    crita <- sqrt(sum(c(newrho - c(rho, sigma), newaa - aa)^2 / (len + 1)^2))
    print(crita)
    if(crita < tol){
        break
    }
}

state <- strtrim(uniloc, 2)
unist <- unique(state)
lunist <- length(unist)
varm <- vector("list")

loclinear <- function(a, i, varm){
    obj <- rep(NA, lunist)
    for(st in 1 : lunist){
        ix <- state == unist[st]
        bb <- log(aggbb[ix, 1])
        Kh <- as.vector(dnorm((AltLong[ix, 1] - AltLong[i, 1])/h[1]) *  dnorm((AltLong[ix, 2] - AltLong[i, 2])/h[2]))
        if(sum(ix) == 1){
            obj[st] <- (bb- a) * Kh * varm[[st]] *( bb- a)
            }else{
        obj[st] <- t(bb - a) %*% diag(Kh) %*% varm[[st]] %*%(bb - a)
        }
    }
    sum(obj)
}
aa <- rep(NA, len)
getaa <- function(i, rho, sigma, v){
    for(st in 1 : length(unique(state))){
    ix <- state == unist[st]
    bb <- sqrt(aggbb[ix, 1])
    
    locdis <- dis[ix, ix]
    mrho <- Matern(locdis, scale = rho, range = sigma, smoothness = v)# 1/(gamma(v) * 2 ^( v- 1)) * (sqrt(2*v) * locdis/rho)^v * besselK(sqrt(2 * v) * locdis/rho, v)
                                        # (1 + sqrt(3) * locdis /rho ) * exp(-sqrt(3) * locdis /rho)#exp(-sqrt(5)  * locdis /  rho) * (1 + sqrt(5) * locdis /rho + 5 * (locdis) ^2 /(3 * rho ^2))
    temp <- (aggbb[ix, 2] * 1/bb)
    Vpara <- temp %*% t(temp) * mrho 
    varm[[st]] <- ginv(Vpara)
    
}
aa<- optim(1, loclinear, gr= NULL, i, method = "L-BFGS-B", lower = -10, upper = 10, control = list(), hessian = FALSE, varm)$par
}
mlike <- function(mrho, aa, v){
    rho <- mrho[1]
    sigma <- mrho[2]
    #v <- mrho[3]
    varm <- Vpara <- vector("list")
    res <- rep(NA, lunist)
    for(st in 1 : length(unique(state))){
    ix <- state == unist[st]
    lix <- sum(ix)
    bb <- sqrt(aggbb[ix, 1])
    
    locdis <- dis[ix, ix]
    mrho <- Matern(locdis, scale =rho, range = sigma, smoothness = v)#1/(gamma(v) * 2 ^( v- 1)) * (sqrt(2*v) * locdis/rho)^v * besselK(sqrt(2 * v) * locdis/rho, v)#(1 + sqrt(3) * locdis /rho ) * exp(-sqrt(3) * locdis /rho)#exp(-sqrt(5)  * locdis /  rho) * (1 + sqrt(5) * locdis /rho + 5 * (locdis) ^2 /(3 * rho ^2))
    temp <- aggbb[ix, 2] * 1/bb
    Vpara[[st]] <- temp %*% t(temp) * mrho# * sigma
    varm[[st]] <- ginv(Vpara[[st]])
#    dev <- bb - aa[ix]
    ## if(lix != 1){
    ## vec1 <- rep(1, lix)
    ## dev <- dev %*% t(dev) - Vpara[[st]]
    ## res[st] <- t(vec1) %*% t(dev) %*% dev %*% vec1
    ## }else{
    ##     dev <- dev^2- Vpara[[st]]
    ##     res[st] <- dev^2
    ## }
}
bm <- bdiag(varm)
   vm <- bdiag(Vpara)
    dev <- log( aggbb[, 1]) - aa
   
    as.numeric((t(dev) %*% bm %*% dev )+ determinant(vm, logarithm = TRUE)[[1]])
#    sum(res)
}

for(j in 1: nsim){
    rho <- newrho[1]
    sigma <- newrho[2]
    v <- 1#newrho[2]
    
    aa <- newaa
    newaa <- do.call(rbind, mclapply(1:len, getaa,rho,  sigma, v, mc.cores = 20))
    newrho <- spg(c(rho, sigma), mlike,gr= NULL,method = 1, project= NULL, lower = c(0.01, 0.01), upper =c(100, 100), projectArgs = NULL, control = list(), quiet = FALSE,   newaa, v )$par
    crita <- sqrt(sum(c(newrho - c(rho, sigma), newaa - aa)^2 / (len + 1)^2))
    print(crita)
    if(crita < tol){
        break
    }
}
lageff <- matrix(NA, 9, 7)
for(i in 1:9){
lageff[i, ] <- apply(loclageff[state%in%div[[i]], ], 2, mean)
}
pdf("lageffave.pdf")
lgtm = 0:6
plot(lageff[1, ] ~ lgtm, col = 1, type = "l", ylim = c(-0.03, 0.05), lty = 1)
for(i in 2:9){
    lines(lageff[i, ] ~ lgtm, col = i, xlab = "Lag time", ylab = "Average lag effect", lty = i)
    }
legend("topright", c("Div1", "Div2", "Div3", "Div4", "Div5", "Div6", "Div7", "Div8", "Div9"), col = 1:9, lty = 1:9)

dev.off()
