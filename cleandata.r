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
load("admissionsHWcontrolDays.RData")
hwind <- dataHW$county.date %in% dataHWcontrol$county.date[dataHWcontrol$hw3]
ctrind <- dataHW$county.date %in% dataHWcontrol$county.date[dataHWcontrol$ctrl3]
dataHW<- cbind(dataHW, hwind, ctrind)
dataHW <- cbind(dataHW, matrix(NA, nrow(dataHW), 7))
subdataHW <- dataHW[dataHW$FIPS == "01001", ]
lag7 <- function(i, loc, subdataHW){
    subdataHW[i, 8:14]<- seq(subdataHW[, 2][i ]-6, subdataHW[, 2][i], 1)%in%  dataHWcontrol$date[dataHWcontrol$hw3 &dataHWcontrol$FIPS == loc]
    subdataHW[i, 8:14]
    
                          
}
temp <- sapply(1 : nrow(subdataHW), lag7, "01001",subdataHW)
subdataHW[, 8:14] <- matrix(unlist(temp), ncol = 7, byrow = T)
subdataHWcontrol <- subdataHW[apply(subdataHW[, 8:14], 1, sum)|subdataHW$ctrind, ]
colnames(subdataHWcontrol)[5] <- "admit.55"
weight <- pmax(3/4 *  (1 - ((0:6)/12)^2), 0)
weight <- weight / sum(weight)
cov1 <- as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight)
cov2 <- log(as.numeric(subdataHWcontrol$date)/365)
Glmres <- glm(admit.55~cov1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
#glmres <- glm(admit.55~ lag1 + cov2 + offset(log(denom)), family = poisson(), data = subdataHWcontrol)
object <- function(bbs){
    sigma <- bbs[4]
    weight <- pmax(3/4 *  (1 - ((0:6)/sigma)^2), 0)
    weight <- weight / sum(weight)
    mean <- exp(bbs[1] + bbs[2] * as.matrix(subdataHWcontrol[, 8:14]) %*% as.matrix(weight) + bbs[3] * cov2 + log(subdataHWcontrol$denom))
    -sum(log(dpois(subdataHWcontrol$admit.55, mean)))
    
}
optim(c(coef(Glmres), 18), object, gr = NULL, method = "L-BFGS-B", lower = c(-10, -2, -10, 11), upper = c(1, 1, 1,  20))
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

getdata <- function(i){
    subdataHW <- dataHW[dataHW$FIPS == uniloc[i], ]
    subdataHW[, 8:14] <- matrix(unlist(sapply(1 : nrow(subdataHW), lag7, uniloc[i], (subdataHW))), ncol = 7, byrow = T)
    subdataHWcontrol <- subdataHW[apply(subdataHW[, 8:14], 1, sum)|subdataHW$ctrind, ]
   aggdata <- data.frame(c(uniloc[i], apply(subdataHWcontrol[, 8:14], 2, mean), (mean(subdataHWcontrol[, 5] /subdataHWcontrol[, 4]))))
    
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
AltLongagg$colorBuckets <- as.numeric(cut(AltLongagg$admit.55, quantile(AltLongagg$admit.55, c(0, seq(0.2, 0.8, length.out = 5),1), type = 3 )))
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