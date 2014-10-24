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
dataHW <-  cbind(data.lag[, 1:4], data.lag$admit.55)
for(i in 2 :7){
load(paste("laggedData/lag", i, ".RData", sep = ""))
dataHW <- rbind(dataHW, cbind(data.lag[, 1:4], data.lag$admit.55))
rm(data.lag)

}
hwind <- dataHW$county.date %in% dataHWcontrol$county.date[dataHWcontrol$hw3]
ctrind <- dataHW$county.date %in% dataHWcontrol$county.date[dataHWcontrol$ctrl3]
dataHW<- cbind(dataHW, hwind, ctrind)
dataHW <- cbind(dataHW, matrix(NA, nrow(dataHW), 7))
subdataHW <- dataHW[dataHW$FIPS == "01001", ]
lag7 <- function(i, loc, subdataHW){
    subdataHW[i, 8:14]<- seq(subdataHW[, 2][i ]-6, subdataHW[, 2][i], 1)%in%  dataHWcontrol$date[dataHWcontrol$hw3 &dataHWcontrol$FIPS == loc]
    subdataHW[i, 8:14]
    
                          
}
sapply(1 : nrow(subdataHW), lag7)
subdataHWcontrol <- subdataHW[apply(subdataHW[, 8:14], 1, sum)|subdataHW$ctrind, ]
weight <- 3/4 *  (1 - ((0:6)/100)^2)
weight <- weight / sum(weight)
cov1 <- as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight)
cov2 <- as.numeric(subdataHWcontrol$date)/365
glmres <- glm(admit.55~cov1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
glmres <- glm(admit.55~ lag1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
object <- function(bbs){
    sigma <- bbs[5]
    weight <- pmax(3/4 *  (1 - ((0:6)/sigma)^2), 0)
  #  weight <- weight / sum(weight)
    mean <- exp(bbs[1] + bbs[2] * as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight) + bbs[3] * cov2 + log(subdataHWcontrol$denom))
    sum((subdataHWcontrol$admit.55 - mean)^2)
}
optim(c(0, 0, 0, 0, 6), object, gr = NULL, method = "L-BFGS-B", lower = c(-10, -2, -2,-2,  6), upper = c(1, 1, 1, 1, 100))
uniloc <- unique(dataHW$FIPS)
len <- length(uniloc)
colnames(dataHW)[5] <- "admit.55"
loccoef <- function(i){
    subdataHW <- dataHW[dataHW$FIPS == uniloc[i], ]
    subdataHW[, 8:14] <- matrix(unlist(sapply(1 : nrow(subdataHW), lag7, uniloc[i], (subdataHW))), ncol = 7, byrow = T)
    subdataHWcontrol <- subdataHW[apply(subdataHW[, 8:14], 1, sum)|subdataHW$ctrind, ]
    weight <- 3/4 *  (1 - ((0:6)/100)^2)
    weight <- weight / sum(weight)
    cov1 <- as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight)
    cov2 <- as.numeric(subdataHWcontrol$date)/365
    glmres <- glm(admit.55~cov1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
    print(i)
    bb <- coef(glmres)
    
}
res <- sapply(1 : 10, loccoef)
