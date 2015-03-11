library(dlnm)
lalg <- function(i){
    return(as.numeric(haix[[ix[i]]]))
}
aLtLong<- sapply(1 : length(ix), lalg)
aLtLong <- t(aLtLong)
dataAltLong <- cbind(dataHWcontrol, aLtLong)
ccs <- read.csv("ccs2.csv")
library(sqldf)
datacss2 <- read.csv(file = "ABCD_07_13_2014.csv")
load("/Users/fjiang/pj6/admissionsHWcontrolDays.RData")
dataHWcontrol$date <-  as.numeric(as.Date(as.character(dataHWcontrol$date), "%Y-%m-%d"))
datacss2$date <- as.numeric(as.Date(as.character(datacss2$date), "%Y-%m-%d"))
#datacss2 <- read.csv.sql(file = "ABCD_07_13_2014.csv", sql = "select * from file where ccs1 = 2 and ccs1 = -999", header = T)
#Link_ID,FIPS,age_gp,Tot_den,Adjusted_tot_den,temperature_daily_county_mean,Weather_total_site,temperature_daily_county_min,temperature_daily_county_max,ccs1,tot_Emerg_admssion,tot_admssion
load("dataHWandCTRLdays.RData")
load(paste("laggedData/lag", 1, ".RData", sep = ""))
dataHW <-  cbind(data.lag[, 1:4], data.lag$admit)
for(i in 2 :7){
load(paste("laggedData/lag", i, ".RData", sep = ""))
dataHW <- rbind(dataHW, cbind(data.lag[, 1:4], data.lag$admit))
rm(data.lag)

}
load("admissionsHWcontrolDays.RData")
hwind <- dataHW$county.date%in% dataHWcontrol$county.date[dataHWcontrol$hw3]
ctrind <- dataHW$county.date %in% dataHWcontrol$county.date[dataHWcontrol$ctrl3]
dataHW<- cbind(dataHW, hwind, ctrind)
dataHW <- cbind(dataHW, matrix(NA, nrow(dataHW), 7))
subdataHW <- dataHW[dataHW$FIPS == "01001", ]
lag7 <- function(i, loc, subdataHW, datahwd, datactrl){
    subdataHW$date <- as.Date(subdataHW$date)
    subdataHW[i, 6:12]<- seq(subdataHW$date[i]-6, subdataHW$date[i], 1)%in% datahwd
    subdataHW[i, 13] <- subdataHW$date[i]%in% datactrl
    subdataHW[i, 6:13]
    
                          

}
temp <- sapply(1 : nrow(subdataHW), lag7, "01001",subdataHW)
subdataHW[, 8:14] <- matrix(unlist(temp), ncol = 7, byrow = T)
subdataHWcontrol <- subdataHW[apply(subdataHW[, 8:14], 1, sum)|subdataHW$ctrind, ]
colnames(subdataHWcontrol)[5] <- "admit"
weight <- pmax(3/4 *  (1 - ((0:6)/12)^2), 0)
weight <- weight / sum(weight)
cov1 <- as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight)
cov2 <- log(as.numeric(subdataHWcontrol$date)/365)
Glmres <- glm(admit~cov1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
#glmres <- glm(admit~ lag1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
object <- function(bbs){
    sigma <- bbs[4]
    weight <- pmax(3/4 *  (1 - ((0:6)/sigma)^2), 0)
    weight <- weight / sum(weight)
    mean <- exp(bbs[1] + bbs[2] * as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight) + bbs[3] * cov2 + log(subdataHWcontrol$denom))
    -sum(log(dpois(subdataHWcontrol$admit, mean)))
    
}
optim(c(coef(Glmres), 18), object, gr = NULL, method = "L-BFGS-B", lower = c(-10, -2, -10, 11), upper = c(1, 1, 1,  20))
uniloc <- unique(dataHWcontrol$FIPS)
len <- length(uniloc)
colnames(dataHW)[5] <- "admit"
loccoef <- function(i){
    load(paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    colnames(subdataHWcontrol)[5] <- "admit"
    weight <- rho
    weight <- weight / sum(weight)
#    cov1 <- as.matrix(subdataHWcontrol[, 6:12]) %*% as.matrix(weight)
    mQ <- as.matrix((subdataHWcontrol[, seq(12, 6, -1)]), ncol = 7)
    cov1 <- crossbasis(mQ, lag = c(0, 6), argvar = list(fun = "lin", cen = TRUE), arglag = list(fun = "poly", degree = df) )
    cov2 <- log(as.numeric(subdataHWcontrol$date - inidate)/365)
    glmres <- glm(admit~cov1  +cov2 +  offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    glmcoef <- summary(glmres)$coefficients
    C <- solve(t(mQ) %*% mQ) %*% t(mQ) %*% cov1
    
    bb <- sum(C %*% glmcoef[2:(df + 2), 1])
    vbb <- C %*% vcov(glmres)[2:(df + 2), 2:(df + 2)] %*% t(C)
    sdbb <- sqrt(t(rep(1, 7)) %*% vbb %*% (rep(1, 7)))
    
    
    
    cbind(bb, sdbb)
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
res <- sapply(1 : 10, loccoef)
aggregate<- aggregate(dataAltLong[, c(4, 6, 14)], by=list(dataAltLong$FIPS), FUN = sum)
AltLong <- apply(dataAltLong[, 25:26], 1, paste, collapse =  ":")
AltLong <- unique(AltLong)
FIPS <- unique(dataAltLong$FIPS)
AltLongagg <- cbind(aggregate, AltLong, FIPS)

gp <- gvisGeoChart(data =AltLongagg, locationvar="AltLong", colorvar = 'hw3', sizevar = "hw3", options= list(region="US", displayMode="markers", sizeAxis="{minValue: 0,  maxSize: 3}",   resolution="provinces", sizeAxis.minSize = 0))
plot(gp)
cut(AltLongagg$hw3, quantile(AltLongagg$hw3, seq(0, 1, 0.2)))

library(maps)
AltLongagg$colorBuckets <- as.numeric(cut(AltLongagg$admit, quantile(AltLongagg$admit, c(0, seq(0.2, 0.8, length.out = 5),1), type = 3 )))
colorsmatched <- AltLongagg$colorBuckets[match(county.fips$fips, AltLongagg$FIPS)]
colors = c("#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", 
    "#980043")
map("county", col = colors[colorsmatched], fill = TRUE, resolution = 0, 
    lty = 0, projection = "polyconic", bg=grey(0.8))
map("state", col = "white", fill = FALSE, add = TRUE, lty = 1, lwd = 1, 
    projection = "polyconic")
map("usa", fill = FALSE, add = TRUE, lty = 1, col = "black", lwd = 1, 
    projection = "polyconic")
title("Hospitalization number for Fluid and electrolyte disorders (CCS 55)")

leg.txt <- c("<3", "3 - 7", "7-13", "13-22", "22-41", ">41")
legend("topright", leg.txt, horiz = TRUE, fill = colors)


library(maps)

AltLongagg$colorBuckets <- as.numeric(cut(AltLongagg$hw3, quantile(AltLongagg$hw3, c(0, seq(0.2, 0.8, length.out = 5),1), type = 3 )))
colorsmatched <- AltLongagg$colorBuckets[match(county.fips$fips, AltLongagg$FIPS)]
colors = c("#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", 
    "#980043")
map("county", col = colors[colorsmatched], fill = TRUE, resolution = 0, 
    lty = 0, projection = "polyconic", bg=grey(0.8))
map("state", col = "white", fill = FALSE, add = TRUE, lty = 1, lwd = 1, 
    projection = "polyconic")
map("usa", fill = FALSE, add = TRUE, lty = 1, col = "black", lwd = 1, 
    projection = "polyconic")
title("Heat weave day counts")

leg.txt <- c("3-12", "12-15", "15-17", "17-18", "18-20", ">20")
legend("topright", leg.txt, horiz = TRUE, fill = colors)






library(maps)

lcumlageff$colorBuckets <- as.numeric(cut((lcumlageff[[1]]), quantile((lcumlageff[[1]]), c(0, seq(0.2, 0.8, length.out = 5),1), type = 3 )))
colorsmatched <- lcumlageff$colorBuckets[match(state.fips$fips, as.numeric(lcumlageff[[2]]))]
colors = c("#F1EEF6", "#D4B9DA", "#C994C7", "#DF65B0", "#DD1C77", 
    "#980043")
#map("county", col = colors[colorsmatched], fill = TRUE, resolution = 0, 
#    lty = 0, projection = "polyconic", bg=grey(0.8))
map("state",col = colors[colorsmatched],   fill = TRUE,lty = 1, lwd = 1, 
    projection = "polyconic", bg=grey(0.8), resolution = 0)
map("usa", fill = FALSE, add = TRUE, lty = 1, col = "black", lwd = 1, 
    projection = "polyconic")
title("Cummulative lag effect")

leg.txt <- c("0.8-1.1", "1.1-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.7", ">1.7")
legend("topright", leg.txt, horiz = TRUE, fill = colors)



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

meanhwday <- aggregate(aggallData[, 6:12], by=list(aggallData$FIPS), FUN = mean)
sumcount <- aggregate(aggallData[, 5], by=list(aggallData$FIPS), FUN = sum)
sumadjust <- aggregate(aggallData[, 5] / aggallData[, 4], by=list(aggallData$FIPS), FUN = sum)
spsz <- aggregate(aggallData[, 5], by=list(aggallData$FIPS), FUN = length)
rho <- rep(NA, 7)
for(i in 1: 7){
    rho[i] <- cor.test(meanhwday[, i + 1], sumadjust[, 2], method = "kendall")$estimate
    }


aggAll1 <- aggregate(cbind(aggallData[, 5]/aggallData[, 4], aggallData[, 6], aggallData[, 8:14]), by=list(aggallData$FIPS), FUN = mean)


colnames(aggAll)[4:10] <- c("lag6", "lag5", "lag4", "lag3", "lag2", "lag1", "lag0" )

tryloc <- function(i){
    res <- try(loccoef(i))
    
    if(class(res) != "try-error"){
        return(res)
    }else{
        return(rep(NA, 2))
        }
    }
aggbb <- lapply(1 : len, tryloc)#, mc.cores = 10)
aggbb <- do.call(rbind, aggbb)
aggbb <- matrix(as.numeric(aggbb), ncol = 3, byrow = T)
aggpvalue <- mclapply(1 : len, tryloc, mc.cores = 10)
aggpvalue <- do.call(rbind, aggpvalue)


locpred <- function(i){
    load(paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    cov1 <- crossbasis(subdataHWcontrol[, seq(12, 6, -1)], lag = c(0, 6))
    cov2 <- log(as.numeric(subdataHWcontrol$date - inidate)/365)
    glmres <- glm(admit~cov1  + cov2 +  offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    mean((exp(predict(glmres)) - subdataHWcontrol$admit)^2, na.rm = T)
}




locmypred <- function(i){
    load(paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    
    weight <- rho
    weight <- weight / sum(weight)
    cov1 <- as.matrix(subdataHWcontrol[, 6:12]) %*% as.matrix(weight)
    cov2 <- log(as.numeric(subdataHWcontrol$date - inidate)/365)
    glmres <- glm(admit~cov1  +cov2 +  offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    print(i)
mean((exp(predict(glmres)) - subdataHWcontrol$admit)^2, na.rm = T)
}


locgetvar <- function(i){
    load(paste("./spatial/hwlagloc55", uniloc[i], sep = "_" ))
    
    weight <- rho
    weight <- weight / sum(weight)
    cov1 <- crossbasis(subdataHWcontrol[, seq(12, 6, -1)], lag = c(0, 6))
   # cov1 <- as.matrix(subdataHWcontrol[, 6:12]) %*% as.matrix(weight)
    cov2 <- log(as.numeric(subdataHWcontrol$date - inidate)/365)
    glmres <- glm(admit~cov1 + cov2 +  offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    print(i)
    coefficients(summary(glmres))[2, 2]
}
ix <- which(sumcount[, 2] > 60)
uniloc <- uniloc[ix]
len <- length(uniloc)
myprederror <- mclapply(1 : len,  locmypred, mc.preschedule = FALSE, mc.cores = 10)
mypreddlnm <- mclapply(1 : len,  locpred, mc.preschedule = FALSE, mc.cores = 10)
myprederror <- as.numeric(unlist(myprederror))
mypreddlnm <- as.numeric(unlist(mypreddlnm))

myvar <- mclapply(1 : len,  locgetvar, mc.preschedule = FALSE, mc.cores = 10)
myvar <- as.numeric(unlist(myvar))



AltLong <- AltLong[ix, ]
AltLong <- as.matrix(AltLong)
h <- c(bw.nrd0(AltLong[, 1]),bw.nrd0(AltLong[, 2])) * len^(-1/15)
dis <- dis[ix, ix]
rho <- aa[4]
myvar <- myvar[ix]
aggbb <- aggbb1[ix, ]
object <- function(aa, s, invSigma, Kh){

 
    estm <- aggbb[, 1] - aa[1] - (s - AltLong) %*% aa[2:3]#apply((aa[, 2 : 3]) *  (s - AltLong), 1, sum)
    temp <- (t(estm) %*% Kh %*% invSigma %*% estm)
    temp
    
}

objectrho <- function(rho,   a0){
    estm <- aggbb[, 1] - a0
    
    
    mrho <-  exp(-sqrt(5)  * dis/10  /  rho) * (1 + sqrt(5) * dis/10 /rho + 5 * (dis/10) ^2 /(3 * rho ^2))
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho
    invSigma <- ginv(Vpara)
   1/2 *  t(estm)  %*% invSigma %*% estm #+1/2 * log(det(Sigma))
    
}

objectrho1 <- function(rho,   a0){
    
    
    
    mrho <- exp(-sqrt(5)  * dis  /  rho) * (1 + sqrt(5) * (dis) /rho + 5 * (dis) ^2 /(3 * rho ^2))
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho
    estm <- aggbb[, 1] %*% t(aggbb[, 1]) -  a0 %*% t(a0) - Vpara
    max(eigen(t(estm)  %*% estm)$values)#+1/2 * log(det(Sigma))
    
}

newaa <- matrix(1, ncol= 3, nrow= len)
temp <- optim(aa, object, gr= NULL, AltLong[i, ], method = "L-BFGS-B", lower = c(-Inf, -Inf, -Inf, 0.1), upper = c(Inf, Inf, Inf, 0.9))
AltLong <- cbind(rep(AltLong[, 1], each = 7), rep(AltLong[, 2], each = 7))
getaa <- function(i, aa, rho, sigma){
 #   print(i)
    s <- AltLong[i, ]
    
    s <- matrix(rep(s, len), ncol = 2, byrow = T)
    Kh <- diag((dnorm((s[, 1] - AltLong[, 1])/ h[1]) * dnorm((s[, 2] - AltLong[, 2])/ h[2])))
    mrho <- exp(-sqrt(5)  * dis /  rho) * (1 + sqrt(5) * dis/rho + 5 * (dis) ^2 /(3 * rho ^2))
    Vpara <- aggbb[, 2] %*% t(aggbb[, 2]) * mrho
    invSigma <- ginv(Vpara)
    aa <- aa[i, ]
    aa <- spg(aa, object, gr = NULL, method = 3, project = NULL, lower = -Inf, upper = Inf, projectArgs = NULL, control = list(trace = TRUE), quiet = FALSE, s, invSigma, Kh )$par
}
mdis <- matrix(NA, 7 * len, 7 * len)
for(l in 1: len){
    for(l1 in 1 : len){
    mdis[ (7 * (l-1) + 1) : (7 * l), (7 * (l1-1) + 1) : (7 * l1) ]  <- dis[l, l1]
    }
    }
estm <- lapply(1:len, getaa, aa, rho1) 
for(j in 1: nsim){
    rho <- newrho
    #sigma <- newrho[2]
    
    aa <- newaa
    newaa <- do.call(rbind, mclapply(1:len, getaa, aa, rho, sigma, mc.cores = 5))
    newrho <- spg(newrho, objectrho1, gr = NULL, method = 3, project = NULL, lower = c(0.05), upper = c( 100), projectArgs = NULL, control = list(trace = FALSE), quiet = FALSE, newaa[, 1])$par
    crita <- sqrt(sum(c(newrho - rho, newaa - aa)^2 / (len + 1)^2))
    print(crita)
    if(crita < tol){
        break
    }
}
tol <- 1e-5
