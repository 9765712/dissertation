## Install packages for the tutorial
#install.packages(c("demography","StMoMo","rgl","googleVis","fanplot", "gdata"))

## Load required libraries
library(demography)
library(StMoMo)
library(rgl)
library(googleVis)
library(fanplot)
library(gdata)
library(MortalityLaws)
library(ggplot2)
library(xlsx)

## Source functions
source("temp_functions.R")

# get data
forecastTime <- 120
ages.fit <- 45:80

cntr = 'RUS'
years <- c(1969:2014)
user = "9765712@gmail.com"
pass = "1563813129"


# load data
russian <- hmd.mx("RUS", user, pass, "Russia")
years <- russian$year
ages<- russian$age
Dxt <- russian$rate[[3]] * russian$pop[[3]]
E0xt <- russian$pop[[3]] + 0.5 * Dxt
Ecxt <- russian$pop[[3]]
qxt <- log(Dxt/E0xt)
# obtain male data
Dxtm <- russian$rate$male * russian$pop$male
E0xtm <- russian$pop$male + 0.5 * Dxtm
Ecxtm <- russian$pop$male
qxtm <- log(Dxtm/E0xtm)

# plot male data

persp3d(ages[0:90], years, qxtm[0:90,], col="skyblue", shade=TRUE,xlab="Ages (0-90)",
        ylab="Years",zlab="Mortality rate (log)")


# obtain female data

Dxtf <- russian$rate$female * russian$pop$female
E0xtf <- russian$pop$female + 0.5 * Dxtf
Ecxtf <- russian$pop$female
qxtf <- log(Dxtf/E0xtf)
# plot female data

persp3d(ages[0:90], years, qxtf[0:90,], col="skyblue", shade=TRUE,xlab="Ages (0-90)",
        ylab="Years",zlab="Mortality rate (log)")

# save data

LUL <- t(qxt[0:100,])
write.xlsx(LUL, "C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/datalul.xlsx")

nbcol = 100
color = rev(rainbow(nbcol, start = 0/6, end = 4/6))
zcol  = cut(qxt[0:100,], nbcol)
persp3d(ages[0:100], years, qxt[0:100,], col=color[zcol], shade=TRUE,xlab="Ages (0-100)",
        ylab="Years",zlab="Mortality rate (log)")

weights <- genWeightMat(ages.fit, years, 3)

## Modelling

modelFit <- array(data = NA, c(8, 3))
colnames(modelFit) <- c("AIC", "BIC", "log likelihood")
rownames(modelFit) <- c("LC", "RH", "CBD", "APC", "M6", "M7", "M8", "PLAT")

## LC model
LC <- lc(link = "log")
LCfit <- fit(LC, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
## plot Lee Carter model fit
plot(LCfit, nCol = 3)

modelFit[1, 1] <- AIC(LCfit)
modelFit[1, 2] <- BIC(LCfit)
modelFit[1, 3] <- LCfit$loglik

# Forecast
LCfor <- forecast(LCfit, h = forecastTime)
LCqxt <- cbind(LCfor$fitted, LCfor$rates)
## get residual fit
LCres <- residuals(LCfit)


## Renshaw and Haberman

RH <- rh(link = "log", cohortAgeFun = "1")
RHfit <- fit(RH, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[2, 1] <- AIC(RHfit)
modelFit[2, 2] <- BIC(RHfit)
modelFit[2, 3] <- RHfit$loglik

RHfor <- forecast(RHfit, h = forecastTime, gc.order = c(1, 1, 0))
RHqxt <- cbind(RHfor$fitted, RHfor$rates)
RHres <- residuals(RHfit)

## APC Currie (2006) - M3

APC <- apc("log")
APCfit <- fit(APC, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[3, 1] <- AIC(APCfit)
modelFit[3, 2] <- BIC(APCfit)
modelFit[3, 3] <- APCfit$loglik

APCfor <- forecast(APCfit, h = forecastTime, gc.order = c(1, 1, 0))
APCqxt <- cbind(APCfor$fitted, APCfor$rates)
APCres <- residuals(APCfit)

## CBD model Cairns (2009) under a Binomial distribution of deaths Haberman and Renshaw (2011) - M5
CBD <- cbd()
CBDfit <- fit(CBD, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[4, 1] <- AIC(CBDfit)
modelFit[4, 2] <- BIC(CBDfit)
modelFit[4, 3] <- CBDfit$loglik

CBDfor <- forecast(CBDfit, h = forecastTime)
CBDqxt <- cbind(CBDfor$fitted, CBDfor$rates)
CBDres <- residuals(CBDfit)

## M6
M6 <- m6()
M6fit <- fit(M6, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[5, 1] <- AIC(M6fit)
modelFit[5, 2] <- BIC(M6fit)
modelFit[5, 3] <- M6fit$loglik

M6for <- forecast(M6fit, h = forecastTime, gc.order = c(2, 0, 0))
M6qxt <- cbind(M6for$fitted, M6for$rates)
M6res <- residuals(M6fit)

## M7 under Binomial setting
M7 <- m7(link = "logit")
M7fit <- fit(M7, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[6, 1] <- AIC(M7fit)
modelFit[6, 2] <- BIC(M7fit)
modelFit[6, 3] <- M7fit$loglik

M7for <- forecast(M7fit, h = forecastTime, gc.order = c(2, 0, 0))
M7qxt <- cbind(M7for$fitted, M7for$rates)
M7res <- residuals(M7fit)

## M8
M8 <- m8(link = "logit", xc = 65)
M8fit <- fit(M8, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[7, 1] <- AIC(M8fit)
modelFit[7, 2] <- BIC(M8fit)
modelFit[7, 3] <- M8fit$loglik

M8for <- forecast(M8fit, h = forecastTime, gc.order = c(2, 0, 0))
M8qxt <- cbind(M8for$fitted, M8for$rates)
M8res <- residuals(M8fit)

## PLAT
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  #\sum g(c)=0, \sum cg(c)=0, \sum c^2g(c)=0
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  #\sum kt[i, ] = 0
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1",
               constFun = constPlat)
PLATfit <- fit(PLAT, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[8, 1] <- AIC(PLATfit)
modelFit[8, 2] <- BIC(PLATfit)
modelFit[8, 3] <- PLATfit$loglik

PLATfor <- forecast(PLATfit, h = forecastTime, gc.order = c(2, 0, 0))
PLATqxt <- cbind(PLATfor$fitted, PLATfor$rates)
PLATres <- residuals(PLATfit)
## Collect fitted models for simulation purposes
modelsFitted <- list(LC = LCfit, RH = RHfit, APC = APCfit, 
               CBD = CBDfit,  M6 = M6fit, M7 = M7fit, M8 = M8fit, PLAT = PLATfit)

## Model fit criteria AIC and BIC
modelFit

## Fit parameter values

par(mfrow=c(1,1))
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parLC1.png")
plot(LCfit, nCol = 3)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parRH1.png")
plot(RHfit, parametricbx = FALSE, nCol = 4)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parAPC1.png")
plot(APCfit, parametricbx = FALSE, nCol = 3)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parCBD1.png")
plot(CBDfit, nCol = 4)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM61.png")
plot(M6fit, parametricbx = FALSE, nCol = 3)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM71.png")
plot(M7fit, parametricbx = FALSE, nCol = 4)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM81.png")
plot(M8fit, parametricbx = FALSE, nCol = 3)
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parPLAT1.png")
plot(PLATfit, parametricbx = FALSE, nCol = 4)
dev.off()

## Plot heatmap residuals
par(mfrow=c(1,1))
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resLC1.png")
plot(LCres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRH1.png")
plot(RHres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resAPC1.png")
plot(APCres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resCBD1.png")
plot(CBDres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM61.png")
plot(M6res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM71.png")
plot(M7res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM81.png")
plot(M8res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLAT1.png")
plot(PLATres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()

## Plot scatterplot residuals
par(mfrow=c(1,1))
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resLC2.png")
plot(LCres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRH2.png")
plot(RHres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resAPC2.png")
plot(APCres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resCBD2.png")
plot(CBDres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM62.png")
plot(M6res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM72.png")
plot(M7res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM82.png")
plot(M8res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
png(filename="C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLAT2.png")
plot(PLATres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

## Plotting - inspection of forecasted mortality rates for 65 year old
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+forecastTime))
plot(years,(Dxt/E0xt)["65", ], pch=21, bg='black',
     xlim = range(1955:2150), ylim = range(0.002, 0.05), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts at age 65, all data")
lines(years_chart, LCqxt["65",], col = "red", lwd = 2)
lines(years_chart, RHqxt["65",], col = "blue", lwd = 2)
lines(years_chart, APCqxt["65",], col = "green", lwd = 2)
lines(years_chart, CBDqxt["65",], col = "lightblue", lwd = 2)
lines(years_chart, M6qxt["65",], col = "pink", lwd = 2)
lines(years_chart, M7qxt["65",], col = "yellow", lwd = 2)
lines(years_chart, M8qxt["65",], col = "purple", lwd = 2)
lines(years_chart, PLATqxt["65",], col = "grey", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("LC","RH", "APC", "CBD", "M6", "M7", "M8", "PLAT"), col=c("red", "blue", "green","lightblue",
       "pink", "yellow", "purple", "grey", "orange"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")

# compare residuals when fitting purely male and female data

# # RENSHAW HABERMAN

# MALE

RH <- rh(link = "log", cohortAgeFun = "1")
RHfitmale <- fit(RH, Dxt = Dxtm, Ext= Ecxtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(RHfitmale)
BIC(RHfitmale)
RHfitmale$loglik

RHformale <- forecast(RHfitmale, h = forecastTime, gc.order = c(1, 1, 0))
RHqxtmale <- cbind(RHformale$fitted, RHformale$rates)
RHresmale <- residuals(RHfitmale)

plot(RHfitmale, parametricbx = FALSE, nCol = 2)
plot(RHresmale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(RHresmale, type = "scatter", reslim = c(-3.5, 3.5))

# FEMALE

RHfitfemale <- fit(RH, Dxt = Dxtf, Ext= Ecxtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(RHfitfemale)
BIC(RHfitfemale)
RHfitfemale$loglik

RHforfemale <- forecast(RHfitfemale, h = forecastTime, gc.order = c(1, 1, 0))
RHqxtfemale <- cbind(RHforfemale$fitted, RHforfemale$rates)
RHresfemale <- residuals(RHfitfemale)

plot(RHfitfemale,parametricbx = FALSE, nCol = 2)
plot(RHresfemale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(RHresfemale, type = "scatter", reslim = c(-3.5, 3.5))




# # M7

# MALE

M7fitmale <- fit(M7, Dxt = Dxtm, Ext= E0xtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(M7fitmale)
BIC(M7fitmale)
M7fitmale$loglik

M7formale <- forecast(M7fitmale, h = forecastTime, gc.order = c(1, 1, 0))
M7qxtmale <- cbind(M7formale$fitted, M7formale$rates)
M7resmale <- residuals(M7fitmale)

plot(M7fitmale, parametricbx = FALSE, nCol = 2)
plot(M7resmale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(M7resmale, type = "scatter", reslim = c(-3.5, 3.5))

# FEMALE

M7fitfemale <- fit(M7, Dxt = Dxtf, Ext= E0xtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(M7fitfemale)
BIC(M7fitfemale)
M7fitfemale$loglik

M7forfemale <- forecast(M7fitfemale, h = forecastTime, gc.order = c(1, 1, 0))
M7qxtfemale <- cbind(M7forfemale$fitted, M7forfemale$rates)
M7resfemale <- residuals(M7fitfemale)

plot(M7fitfemale,parametricbx = FALSE, nCol = 2)
plot(M7resfemale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(M7resfemale, type = "scatter", reslim = c(-3.5, 3.5))


# # PLAT

# MALE

PLATfitmale <- fit(PLAT, Dxt = Dxtm, Ext= E0xtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

PLATformale <- forecast(PLATfitmale, h = forecastTime, gc.order = c(2, 0, 0))
PLATqxtmale <- cbind(PLATformale$fitted, PLATformale$rates)
PLATresmale <- residuals(PLATfitmale)

plot(PLATfitmale, parametricbx = FALSE,  nCol = 2)
plot(PLATresmale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(PLATresmale, type = "scatter", reslim = c(-3.5, 3.5))

# FEMALE

PLATfitfemale <- fit(PLAT, Dxt = Dxtf, Ext= E0xtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

PLATforfemale <- forecast(PLATfitfemale, h = forecastTime, gc.order = c(2, 0, 0))
PLATqxtfemale <- cbind(PLATforfemale$fitted, PLATforfemale$rates)
PLATresfemale <- residuals(PLATfitfemale)

plot(PLATfitfemale, parametricbx = FALSE, nCol = 2)
plot(PLATresfemale, type = "colourmap", reslim = c(-3.5, 3.5))
plot(PLATresfemale, type = "scatter", reslim = c(-3.5, 3.5))


# plot forecasts by gender

# males

years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+forecastTime))
plot(years,(Dxtm/E0xtm)["65", ], pch=21, bg='black',
     xlim = range(1955:2140), ylim = range(0.005, 0.065), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts at age 65, male data only")
lines(years_chart, RHqxtmale["65",], col = "blue", lwd = 2)
lines(years_chart, M7qxtmale["65",], col = "red", lwd = 2)
lines(years_chart, PLATqxtmale["65",], col = "green", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("RH", "M7", "PLAT"), col=c("blue", "red", "green"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")


# females

years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+forecastTime))
plot(years,(Dxtf/E0xtf)["65", ], pch=21, bg='black',
     xlim = range(1955:2140), ylim = range(0.005, 0.028), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts at age 65, female data only")
lines(years_chart, RHqxtfemale["65",], col = "blue", lwd = 2)
lines(years_chart, M7qxtfemale["65",], col = "red", lwd = 2)
lines(years_chart, PLATqxtfemale["65",], col = "green", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("RH", "M7", "PLAT"), col=c("blue", "red", "green"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")






## Extrapolate data to omega age = 120
LCextrapolate <- kannistoExtrapolation(LCqxt, ages.fit, years_chart)
LCqxtExtr <- LCextrapolate$qxt

RHextrapolate <- kannistoExtrapolation(RHqxt, ages.fit, years_chart)
RHqxtExtr <- RHextrapolate$qxt

APCextrapolate <- kannistoExtrapolation(APCqxt, ages.fit, years_chart)
APCqxtExtr <- APCextrapolate$qxt

CBDextrapolate <- kannistoExtrapolation(CBDqxt, ages.fit, years_chart)
CBDqxtExtr <- CBDextrapolate$qxt

M6extrapolate <- kannistoExtrapolation(M6qxt, ages.fit, years_chart)
M6qxtExtr <- M6extrapolate$qxt

M7extrapolate <- kannistoExtrapolation(M7qxt, ages.fit, years_chart)
M7qxtExtr <- M7extrapolate$qxt

M8extrapolate <- kannistoExtrapolation(M8qxt, ages.fit, years_chart)
M8qxtExtr <- M8extrapolate$qxt

extrapolate <- kannistoExtrapolation(PLATqxt, ages.fit, years_chart)
PLATqxtExtr <- extrapolate$qxt

models <- list(LCqxtExtr = LCqxtExtr, RHqxtExtr = RHqxtExtr, APCqxtExtr = APCqxtExtr, 
               CBDqxtExtr = CBDqxtExtr,  M6qxtExtr = M6qxtExtr, 
               M7qxtExtr = M7qxtExtr, M8qxtExtr = M8qxtExtr, PLATqxtExtr = PLATqxtExtr)

## Annuity projection
## Assumption annuity due deffered at 65

## Read in portfolio data
portfolio <- read.csv('portfolio_2.csv')
# portfolio <- read.csv('portfolio.csv')
portfolio
## Read in experience data
experience.factors <- read.csv('experience-factors.csv')

ages.fit <- 25:120
valyear <- 2015
pension <- 120000

# calculate ages of insured
portfolio$age <- valyear - portfolio$YoB

experience.factors$total <- (experience.factors$Male + experience.factors$Female)/2

expF <- experience.factors$total[ages.fit]

BEL <- array(NA, c(7,1))
rownames(BEL) <- c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT")
colnames(BEL) <- "BEL"

for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
    output2[[i]] <- DFcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
  }
  BEL[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}


## Show BEL per model
df = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL)
Column <- gvisColumnChart(df, options=list(series="[{color:'#0A79BF'}]"))
plot(Column)
Column
## SIMULATIONS for SCR calculation purposes
set.seed(1234)

# number of simulations
nsim <- 200
models2run <- length(modelsFitted)
ages.fit <- 25:90

modelSim <- list()
selectBEL <- list()

for (m in 1:models2run){
  modelSim[[m]] <- simulate(modelsFitted[[m]], nsim = nsim, h = forecastTime)
  collectBEL <- array(NA, c(200,1))
  for (s in 1:nsim){
    prem_s <- 0
    ben_s <- 0
    qx <- cbind(modelSim[[m]]$fitted[, , s], modelSim[[m]]$rates[, , s])
    extrapolate <- kannistoExtrapolation(qx, ages.fit, years_chart)
    for (i in 1:nrow(portfolio)){
      ben_s <- ben_s + DFcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                  pensionAge = 65, valyear = valyear, ir = 0.02, type = 1)*pension
      prem_s <- prem_s + DFcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                    pensionAge = 65, valyear = valyear, ir = 0.02, type = 2)*portfolio$Premium[i]
    }
    collectBEL[s, 1] <- ben_s - prem_s
  }
  selectBEL[[m]] <- quantile(collectBEL, probs = 0.995, type = 1)
}

SCR <- as.numeric(selectBEL) - BEL
colnames(SCR) <- "SCR"

## Plot simulations for LC model
qxt <- Dxt/E0xt
plot(LCfit$years, qxt["60", ], xlim = c(1950, 2129), ylim = range(0.0048, 0.1),
     xlab = "Years", ylab = "Mortality rates", main = "Mortality rates at age 60",
     pch = 20, log = "y", type = "l")
matlines(modelSim[[1]]$years, modelSim[[1]]$rates["60", , 1:20], type = "l", lty = 1, col = 1:20)

## Plot model uncertainity
probs <- c(0.5, 2.5, 10, 25, 75, 90, 97.5, 99.5)
plot(LCfit$years, qxt["60", ], xlim = c(1950, 2129), ylim = c(0.001, 0.1),
     xlab = "Years", ylab = "Mortality rates at age 60", main = "Uncertainity associated with a model forecast",
     pch = 20, log = "y")
fan(t(modelSim[[1]]$rates["60", , ]), start = 2010, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("yellow", "darkgreen")), ln = NULL)

## Show SCR per model
df = data.frame(models = c("LC", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL, SCR = SCR)
Column <- gvisColumnChart(df,
                          options=list(series="[{color:'#0A79BF'}, {color: '#B2246B'}]"))
plot(Column)
