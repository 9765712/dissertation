# FITTING AND GOODNESS OF FIT ANALYSIS

library(demography)
library(StMoMo)
library(MortalityLaws)
library(xlsx)

setwd("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/R/Dissertation")

source("functions.R")

horizon <- 40
ages.fit <- 45:80

user = "9765712@gmail.com"
pass = "1563813129"

# load combined data

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

# obtain female data

Dxtf <- russian$rate$female * russian$pop$female
E0xtf <- russian$pop$female + 0.5 * Dxtf
Ecxtf <- russian$pop$female
qxtf <- log(Dxtf/E0xtf)

weights <- genWeightMat(ages.fit, years, 3)

# save data for plotting

# data <- t(qxt[0:100,])
# write.xlsx(LUL, "C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/data.xlsx")

# data_m <- t(qxtm[0:100,])
# write.xlsx(data_m, "C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/datam.xlsx")

# data_f <- t(qxtf[0:100,])
# write.xlsx(data_f, "C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/dataf.xlsx")

## Modelling and saving results

modelFit <- array(data = NA, c(8, 3))
colnames(modelFit) <- c("AIC", "BIC", "log likelihood")
rownames(modelFit) <- c("LC", "RH", "APC", "CBD", "M6", "M7", "M8", "PLAT")

## LC model

LC <- lc(link = "log")
LCfit <- fit(LC, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[1, 1] <- AIC(LCfit)
modelFit[1, 2] <- BIC(LCfit)
modelFit[1, 3] <- LCfit$loglik

LCfor <- forecast(LCfit, h = horizon)
LCqxt <- cbind(LCfor$fitted, LCfor$rates)
LCres <- residuals(LCfit)

## Renshaw and Haberman

RH <- rh(link = "log", cohortAgeFun = "1")
RHfit <- fit(RH, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[2, 1] <- AIC(RHfit)
modelFit[2, 2] <- BIC(RHfit)
modelFit[2, 3] <- RHfit$loglik

RHfor <- forecast(RHfit, h = horizon, gc.order = c(1, 1, 0))
RHqxt <- cbind(RHfor$fitted, RHfor$rates)
RHres <- residuals(RHfit)

## APC

APC <- apc("log")
APCfit <- fit(APC, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[3, 1] <- AIC(APCfit)
modelFit[3, 2] <- BIC(APCfit)
modelFit[3, 3] <- APCfit$loglik

APCfor <- forecast(APCfit, h = horizon, gc.order = c(1, 1, 0))
APCqxt <- cbind(APCfor$fitted, APCfor$rates)
APCres <- residuals(APCfit)

## CBD

CBD <- cbd()
CBDfit <- fit(CBD, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[4, 1] <- AIC(CBDfit)
modelFit[4, 2] <- BIC(CBDfit)
modelFit[4, 3] <- CBDfit$loglik

CBDfor <- forecast(CBDfit, h = horizon)
CBDqxt <- cbind(CBDfor$fitted, CBDfor$rates)
CBDres <- residuals(CBDfit)

## M6

M6 <- m6()
M6fit <- fit(M6, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[5, 1] <- AIC(M6fit)
modelFit[5, 2] <- BIC(M6fit)
modelFit[5, 3] <- M6fit$loglik

M6for <- forecast(M6fit, h = horizon, gc.order = c(2, 0, 0))
M6qxt <- cbind(M6for$fitted, M6for$rates)
M6res <- residuals(M6fit)

## M7

M7 <- m7(link = "logit")
M7fit <- fit(M7, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[6, 1] <- AIC(M7fit)
modelFit[6, 2] <- BIC(M7fit)
modelFit[6, 3] <- M7fit$loglik

M7for <- forecast(M7fit, h = horizon, gc.order = c(2, 0, 0))
M7qxt <- cbind(M7for$fitted, M7for$rates)
M7res <- residuals(M7fit)

## M8

M8 <- m8(link = "logit", xc = 65)
M8fit <- fit(M8, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[7, 1] <- AIC(M8fit)
modelFit[7, 2] <- BIC(M8fit)
modelFit[7, 3] <- M8fit$loglik

M8for <- forecast(M8fit, h = horizon, gc.order = c(2, 0, 0))
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
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
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

PLATfor <- forecast(PLATfit, h = horizon, gc.order = c(2, 0, 0))
PLATqxt <- cbind(PLATfor$fitted, PLATfor$rates)
PLATres <- residuals(PLATfit)

modelsFitted <- list(LC = LCfit, RH = RHfit, APC = APCfit, 
               CBD = CBDfit,  M6 = M6fit, M7 = M7fit, M8 = M8fit, PLAT = PLATfit)

modelFit

## Plot parameters

par(mfrow=c(1,1))
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parLC1.pdf", width=7, height=2.7)
plot(LCfit, nCol = 3)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parRH1.pdf", width=7, height=2.7)
plot(RHfit, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parAPC1.pdf", width=7, height=2.7)
plot(APCfit, nCol = 3)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parCBD1.pdf", width=7, height=2.7)
plot(CBDfit, parametricbx = FALSE, nCol = 2)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM61.pdf", width=7, height=2.7)
plot(M6fit, parametricbx = FALSE, nCol = 3)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM71.pdf", width=7, height=2.7)
plot(M7fit, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM81.pdf", width=7, height=2.7)
plot(M8fit, parametricbx = FALSE, nCol = 3)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parPLAT1.pdf", width=7, height=2.7)
plot(PLATfit, parametricbx = FALSE, nCol = 4)
dev.off()

## Plot heatmap residuals

par(mfrow=c(1,1))
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resLC1.pdf", width=5, height=5)
plot(LCres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRH1.pdf", width=5, height=5)
plot(RHres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resAPC1.pdf", width=5, height=5)
plot(APCres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resCBD1.pdf", width=5, height=5)
plot(CBDres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM61.pdf", width=5, height=5)
plot(M6res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM71.pdf", width=5, height=5)
plot(M7res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM81.pdf", width=5, height=5)
plot(M8res, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLAT1.pdf", width=5, height=5)
plot(PLATres, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()

## Plot scatterplot residuals

par(mfrow=c(1,1))
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resLC2.pdf", width=7, height=2.7)
plot(LCres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRH2.pdf", width=7, height=2.7)
plot(RHres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resAPC2.pdf", width=7, height=2.7)
plot(APCres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resCBD2.pdf", width=7, height=2.7)
plot(CBDres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM62.pdf", width=7, height=2.7)
plot(M6res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM72.pdf", width=7, height=2.7)
plot(M7res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM82.pdf", width=7, height=2.7)
plot(M8res, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLAT2.pdf", width=7, height=2.7)
plot(PLATres, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

## Forecast of mortality for age 65

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/forecast65.pdf", width=7, height=5)
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+horizon))
plot(years,(Dxt/E0xt)["65", ], pch=21, bg='black',
     xlim = range(1959:2053), ylim = range(0.002, 0.045), xlab = "Years", ylab = "Mortality rate at age 65",
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
dev.off()

## Compare models for male and female data separately

## Renshaw-Haberman

# MALE

RH <- rh(link = "log", cohortAgeFun = "1")
RHfitmale <- fit(RH, Dxt = Dxtm, Ext= Ecxtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(RHfitmale)
BIC(RHfitmale)
RHfitmale$loglik

RHformale <- forecast(RHfitmale, h = horizon, gc.order = c(1, 1, 0))
RHqxtmale <- cbind(RHformale$fitted, RHformale$rates)
RHresmale <- residuals(RHfitmale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parRHmale.pdf", width=7, height=2.7)
plot(RHfitmale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRHmale1.pdf", width=5, height=5)
plot(RHresmale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRHmale2.pdf", width=7, height=2.7)
plot(RHresmale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

# FEMALE

RHfitfemale <- fit(RH, Dxt = Dxtf, Ext= Ecxtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(RHfitfemale)
BIC(RHfitfemale)
RHfitfemale$loglik

RHforfemale <- forecast(RHfitfemale, h = horizon, gc.order = c(1, 1, 0))
RHqxtfemale <- cbind(RHforfemale$fitted, RHforfemale$rates)
RHresfemale <- residuals(RHfitfemale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parRHfemale.pdf", width=7, height=2.7)
plot(RHfitfemale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRHfemale1.pdf", width=5, height=5)
plot(RHresfemale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resRHfemale2.pdf", width=7, height=2.7)
plot(RHresfemale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

## M7

# MALE

M7fitmale <- fit(M7, Dxt = Dxtm, Ext= E0xtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(M7fitmale)
BIC(M7fitmale)
M7fitmale$loglik

M7formale <- forecast(M7fitmale, h = horizon, gc.order = c(1, 1, 0))
M7qxtmale <- cbind(M7formale$fitted, M7formale$rates)
M7resmale <- residuals(M7fitmale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM7male.pdf", width=7, height=2.7)
plot(M7fitmale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM7male1.pdf", width=5, height=5)
plot(M7resmale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM7male2.pdf", width=7, height=2.7)
plot(M7resmale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

# FEMALE

M7fitfemale <- fit(M7, Dxt = Dxtf, Ext= E0xtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

AIC(M7fitfemale)
BIC(M7fitfemale)
M7fitfemale$loglik

M7forfemale <- forecast(M7fitfemale, h = horizon, gc.order = c(1, 1, 0))
M7qxtfemale <- cbind(M7forfemale$fitted, M7forfemale$rates)
M7resfemale <- residuals(M7fitfemale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parM7female.pdf", width=7, height=2.7)
plot(M7fitfemale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM7female1.pdf", width=5, height=5)
plot(M7resfemale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resM7female2.pdf", width=7, height=2.7)
plot(M7resfemale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

## PLAT

# MALE

PLATfitmale <- fit(PLAT, Dxt = Dxtm, Ext= E0xtm, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

PLATformale <- forecast(PLATfitmale, h = horizon, gc.order = c(2, 0, 0))
PLATqxtmale <- cbind(PLATformale$fitted, PLATformale$rates)
PLATresmale <- residuals(PLATfitmale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parPLATmale.pdf", width=7, height=2.7)
plot(PLATfitmale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLATmale1.pdf", width=5, height=5)
plot(PLATresmale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLATmale2.pdf", width=7, height=2.7)
plot(PLATresmale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

# FEMALE

PLATfitfemale <- fit(PLAT, Dxt = Dxtf, Ext= E0xtf, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

PLATforfemale <- forecast(PLATfitfemale, h = horizon, gc.order = c(2, 0, 0))
PLATqxtfemale <- cbind(PLATforfemale$fitted, PLATforfemale$rates)
PLATresfemale <- residuals(PLATfitfemale)

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/parPLATfemale.pdf", width=7, height=2.7)
plot(PLATfitfemale, parametricbx = FALSE, nCol = 4)
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLATfemale1.pdf", width=5, height=5)
plot(PLATresfemale, type = "colourmap", reslim = c(-3.5, 3.5))
dev.off()
pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/resPLATfemale2.pdf", width=7, height=2.7)
plot(PLATresfemale, type = "scatter", reslim = c(-3.5, 3.5))
dev.off()

## Plot forecasts for each gender

# Males

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/forecast65male.pdf", width=7, height=5)
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+horizon))
plot(years,(Dxtm/E0xtm)["65", ], pch=21, bg='black',
     xlim = range(1959:2053), ylim = range(0.005, 0.065), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts at age 65, male data only")
lines(years_chart, RHqxtmale["65",], col = "blue", lwd = 2)
lines(years_chart, M7qxtmale["65",], col = "yellow", lwd = 2)
lines(years_chart, PLATqxtmale["65",], col = "gray", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("RH", "M7", "PLAT"), col=c("blue", "yellow", "gray"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")
dev.off()

# Females

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/forecast65female.pdf", width=7, height=5)
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+horizon))
plot(years,(Dxtf/E0xtf)["65", ], pch=21, bg='black',
     xlim = range(1959:2053), ylim = range(0.005, 0.025), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts at age 65, female data only")
lines(years_chart, RHqxtfemale["65",], col = "blue", lwd = 2)
lines(years_chart, M7qxtfemale["65",], col = "yellow", lwd = 2)
lines(years_chart, PLATqxtfemale["65",], col = "gray", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("RH", "M7", "PLAT"), col=c("blue", "yellow", "gray"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")
dev.off()

## Robustness analysis

# extract data for years of market economy, 1992 - 2014

russian_mkt <- extract.years(russian, 1992:2014)

# obtain recent data

years_m <- russian_mkt$year
ages_m<- russian_mkt$age
Dxt_m <- russian_mkt$rate[[3]] * russian_mkt$pop[[3]]
E0xt_m <- russian_mkt$pop[[3]] + 0.5 * Dxt_m
Ecxt_m <- russian_mkt$pop[[3]]
qxt_m <- log(Dxt_m/E0xt_m)

weights_m <- genWeightMat(ages.fit, 1992:2014, 3)
horizon_m <- 20

modelFit_m <- array(data = NA, c(8, 3))
colnames(modelFit_m) <- c("AIC", "BIC", "log likelihood")
rownames(modelFit_m) <- c("LC", "RH", "APC", "CBD", "M6", "M7", "M8", "PLAT")

# LC

LC_m <- lc(link = "log")
LCfit_m <- fit(LC_m, Dxt = Dxt_m, Ext= Ecxt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)

modelFit_m[1, 1] <- AIC(LCfit_m)
modelFit_m[1, 2] <- BIC(LCfit_m)
modelFit_m[1, 3] <- LCfit_m$loglik

LCfor_m <- forecast(LCfit_m, h = horizon_m)
LCqxt_m <- cbind(LCfor_m$fitted, LCfor_m$rates)
LCres_m <- residuals(LCfit_m)

# RH

RH_m <- rh(link = "log", cohortAgeFun = "1")
RHfit_m <- fit(RH_m, Dxt = Dxt_m, Ext= Ecxt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)

modelFit_m[2, 1] <- AIC(RHfit_m)
modelFit_m[2, 2] <- BIC(RHfit_m)
modelFit_m[2, 3] <- RHfit_m$loglik

RHfor_m <- forecast(RHfit_m, h = horizon_m, gc.order = c(1, 1, 0))
RHqxt_m <- cbind(RHfor_m$fitted, RHfor_m$rates)
RHres_m <- residuals(RHfit_m)

plot(RHfit_m, parametricbx = FALSE, nCol = 2)
plot(RHres_m, type = "colourmap", reslim = c(-3.5, 3.5))
plot(RHres_m, type = "scatter", reslim = c(-3.5, 3.5))

# APC

APC_m <- apc("log")
APCfit_m <- fit(APC_m, Dxt = Dxt_m, Ext= Ecxt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m, start.ax = LCfit_m$ax,
              start.bx = LCfit_m$bx, start.kt = LCfit_m$kx)

modelFit_m[3, 1] <- AIC(APCfit_m)
modelFit_m[3, 2] <- BIC(APCfit_m)
modelFit_m[3, 3] <- APCfit_m$loglik

APCfor_m <- forecast(APCfit_m, h = horizon_m, gc.order = c(1, 1, 0))
APCqxt_m <- cbind(APCfor_m$fitted, APCfor_m$rates)
APCres_m <- residuals(APCfit_m)

# CBD

CBD_m <- cbd()
CBDfit_m <- fit(CBD_m, Dxt = Dxt_m, Ext= E0xt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)

modelFit_m[4, 1] <- AIC(CBDfit_m)
modelFit_m[4, 2] <- BIC(CBDfit_m)
modelFit_m[4, 3] <- CBDfit_m$loglik

CBDfor_m <- forecast(CBDfit_m, h = horizon_m)
CBDqxt_m <- cbind(CBDfor_m$fitted, CBDfor_m$rates)
CBDres_m <- residuals(CBDfit_m)

# M6 

M6_m <- m6()
M6fit_m <- fit(M6_m, Dxt = Dxt_m, Ext= E0xt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)
modelFit_m[5, 1] <- AIC(M6fit_m)
modelFit_m[5, 2] <- BIC(M6fit_m)
modelFit_m[5, 3] <- M6fit_m$loglik

M6for_m <- forecast(M6fit_m, h = horizon_m, gc.order = c(2, 0, 0))
M6qxt_m <- cbind(M6for_m$fitted, M6for_m$rates)
M6res_m <- residuals(M6fit_m)

# M7

M7fit_m <- fit(M7, Dxt = Dxt_m, Ext= Ecxt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)

AIC(M7fit_m)
BIC(M7fit_m)
M7fit_m$loglik

M7for_m <- forecast(M7fit_m, h = horizon_m, gc.order = c(1, 1, 0))
M7qxt_m <- cbind(M7for_m$fitted, M7for_m$rates)
M7res_m <- residuals(M7fit_m)

plot(M7fit_m, parametricbx = FALSE, nCol = 2)
plot(M7res_m, type = "colourmap", reslim = c(-3.5, 3.5))
plot(M7res_m, type = "scatter", reslim = c(-3.5, 3.5))

# M8

M8_m <- m8(link = "logit", xc = 65)
M8fit_m <- fit(M8_m, Dxt = Dxt_m, Ext= E0xt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)
modelFit_m[7, 1] <- AIC(M8fit_m)
modelFit_m[7, 2] <- BIC(M8fit_m)
modelFit_m[7, 3] <- M8fit_m$loglik

M8for_m <- forecast(M8fit_m, h = horizon_m, gc.order = c(2, 0, 0))
M8qxt_m <- cbind(M8for_m$fitted, M8for_m$rates)
M8res_m <- residuals(M8fit_m)

# PLAT

PLAT_m <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1",
               constFun = constPlat)
PLATfit_m <- fit(PLAT, Dxt = Dxt_m, Ext= E0xt_m, ages = ages, years = years_m, ages.fit = ages.fit, wxt = weights_m)

modelFit_m[8, 1] <- AIC(PLATfit_m)
modelFit_m[8, 2] <- BIC(PLATfit_m)
modelFit_m[8, 3] <- PLATfit_m$loglik

PLATfor_m <- forecast(PLATfit_m, h = horizon_m, gc.order = c(2, 0, 0))
PLATqxt_m <- cbind(PLATfor_m$fitted, PLATfor_m$rates)
PLATres_m <- residuals(PLATfit_m)

# Plot forecast at 65

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/recentforecast65.pdf", width=7, height=5)
years_chart2 <- c(years_m, (years_m[length(years_m)]+1):(years_m[length(years_m)]+horizon_m))
plot(years_m,(Dxt_m/E0xt_m)["65", ], pch=21, bg='black',
     xlim = range(1992:2034), ylim = range(0.003, 0.04), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Forecasts at age 65, recent data")
lines(years_chart2, LCqxt_m["65",], col = "red", lwd = 2)
lines(years_chart2, RHqxt_m["65",], col = "blue", lwd = 2)
lines(years_chart2, APCqxt_m["65",], col = "green", lwd = 2)
lines(years_chart2, CBDqxt_m["65",], col = "lightblue", lwd = 2)
lines(years_chart2, M6qxt_m["65",], col = "pink", lwd = 2)
lines(years_chart2, M7qxt_m["65",], col = "yellow", lwd = 2)
lines(years_chart2, M8qxt_m["65",], col = "purple", lwd = 2)
lines(years_chart2, PLATqxt_m["65",], col = "grey", lwd = 2)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("LC","RH", "APC", "CBD", "M6", "M7", "M8", "PLAT"), col=c("red", "blue", "green","lightblue",
                                                                                 "pink", "yellow", "purple", "grey", "orange"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")
dev.off()