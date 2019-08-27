# SIMULATION

# devtools::install_github("TARF/insureR")

library(demography)
library(StMoMo)
library(MortalityLaws)

setwd("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/R/Dissertation")

source("functions.R")

# get data

horizon <- 120
ages.fit <- 45:80

cntr = 'RUS'
years <- c(1992:2014)
user = "9765712@gmail.com"
pass = "1563813129"

# load data
russian <- hmd.mx("RUS", user, pass, "Russia")
russian <- extract.years(russian, 1992:2014)
years <- russian$year
ages<- russian$age
Dxt <- russian$rate[[3]] * russian$pop[[3]]
E0xt <- russian$pop[[3]] + 0.5 * Dxt
Ecxt <- russian$pop[[3]]
qxt <- log(Dxt/E0xt)

weights <- genWeightMat(ages.fit, years, 3)

## Modelling

modelFit <- array(data = NA, c(8, 3))
colnames(modelFit) <- c("AIC", "BIC")
rownames(modelFit) <- c("LC", "RH", "APC", "CBD", "M6", "M7", "M8", "PLAT")

## Lee-Carter

LC <- lc(link = "log")
LCfit <- fit(LC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[1, 1] <- AIC(LCfit)
modelFit[1, 2] <- BIC(LCfit)
modelFit[1, 3] <- LCfit$loglik

LCfor <- forecast(LCfit, h = horizon)
LCqxt <- cbind(LCfor$fitted, LCfor$rates)
LCres <- residuals(LCfit)

# Renshaw-Haberman

RH <- rh(link = "log", cohortAgeFun = "1")
RHfit <- fit(RH, Dxt = Dxt, Ext= Ecxt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[2, 1] <- AIC(RHfit)
modelFit[2, 2] <- BIC(RHfit)
modelFit[2, 3] <- RHfit$loglik

RHfor <- forecast(RHfit, h = horizon, gc.order = c(1, 1, 0))
RHqxt <- cbind(RHfor$fitted, RHfor$rates)
RHres <- residuals(RHfit)

# Age-Period-Cohort

APC <- apc("log")
APCfit <- fit(APC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[3, 1] <- AIC(APCfit)
modelFit[3, 2] <- BIC(APCfit)
modelFit[3, 3] <- APCfit$loglik

APCfor <- forecast(APCfit, h = horizon, gc.order = c(1, 1, 0))
APCqxt <- cbind(APCfor$fitted, APCfor$rates)
APCres <- residuals(APCfit)

# Cairns-Blake-Dowd

CBD <- cbd()
CBDfit <- fit(CBD, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[4, 1] <- AIC(CBDfit)
modelFit[4, 2] <- BIC(CBDfit)
modelFit[4, 3] <- CBDfit$loglik

CBDfor <- forecast(CBDfit, h = horizon)
CBDqxt <- cbind(CBDfor$fitted, CBDfor$rates)
CBDres <- residuals(CBDfit)

# M6 (CBD with cohort effect)

M6 <- m6()
M6fit <- fit(M6, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[5, 1] <- AIC(M6fit)
modelFit[5, 2] <- BIC(M6fit)
modelFit[5, 3] <- M6fit$loglik

M6for <- forecast(M6fit, h = horizon, gc.order = c(2, 0, 0))
M6qxt <- cbind(M6for$fitted, M6for$rates)
M6res <- residuals(M6fit)

# M7 (CBD with quadratic age effect)

M7 <- m7(link = "logit")
M7fit <- fit(M7, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[6, 1] <- AIC(M7fit)
modelFit[6, 2] <- BIC(M7fit)
modelFit[6, 3] <- M7fit$loglik

M7for <- forecast(M7fit, h = horizon, gc.order = c(2, 0, 0))
M7qxt <- cbind(M7for$fitted, M7for$rates)
M7res <- residuals(M7fit)

# M8 (CBD with reducing cohort effect)

M8 <- m8(link = "logit", xc = 65)
M8fit <- fit(M8, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[7, 1] <- AIC(M8fit)
modelFit[7, 2] <- BIC(M8fit)
modelFit[7, 3] <- M8fit$loglik

M8for <- forecast(M8fit, h = horizon, gc.order = c(2, 0, 0))
M8qxt <- cbind(M8for$fitted, M8for$rates)
M8res <- residuals(M8fit)

# PLAT

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

# Forecast of mortality rates for age 65

years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+horizon))
plot(years,(Dxt/E0xt)["65", ], pch=21, col='blue', bg='lightblue',
     xlim = range(1990:2150), ylim = range(0.01, 0.045), xlab = "Years projection", ylab = "Mortality rate at age 65",
     main = "Different models' forecasts", bty="n")
lines(years_chart, LCqxt["65",], col = "red", lwd = 2)
lines(years_chart, RHqxt["65",], col = "orange", lwd = 2)
lines(years_chart, APCqxt["65",], col = "green", lwd = 2)
lines(years_chart, CBDqxt["65",], col = "lightblue", lwd = 2)
lines(years_chart, M6qxt["65",], col = "pink", lwd = 2)
lines(years_chart, M7qxt["65",], col = "yellow", lwd = 2)
lines(years_chart, M8qxt["65",], col = "purple", lwd = 2)
lines(years_chart, PLATqxt["65",], col = "grey", lwd = 2)
abline(h = 0, v = 2014, col = "gray60")
legend("topright", c("LC", "RH", "APC", "CBD", "M6", "M7", "M8", "PLAT"), col=c("red", "orange",  "green","lightblue",
       "pink", "yellow", "purple", "grey"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")

# Kannisto-extrapolatiton of data to age 120

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

# Annuity due at 65 y/o

portfolio <- read.csv('sberbank.csv')

ages.fit <- 45:120
valyear <- 2015
pension <- 120000

portfolio$age <- valyear - portfolio$YoB

# Calculate BEL

BEL <- array(NA, c(8,1))
rownames(BEL) <- c("LC","RH", "APC", "CBD", "M6", "M7", "M8", "PLAT")
colnames(BEL) <- "BEL"

for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- DFcashflow(models[[m]], ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.08, type = 1)*pension
    output2[[i]] <- DFcashflow(models[[m]], ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, valyear = valyear, ir = 0.08, type = 2)*portfolio$Premium[i]
  }
  BEL[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}

# Simulations for SCR calculation, 80%, 90% and 99.5%

set.seed(1234)

nsim <- 200
models2run <- length(modelsFitted)
ages.fit <- 45:80

modelSim <- list()

selectBEL80 <-list()
selectBEL90 <- list()
selectBEL995 <- list()

for (m in 1:models2run){
  modelSim[[m]] <- simulate(modelsFitted[[m]], nsim = nsim, h = horizon)
  collectBEL <- array(NA, c(200,1))
  for (s in 1:nsim){
    prem_s <- 0
    ben_s <- 0
    qx <- cbind(modelSim[[m]]$fitted[, , s], modelSim[[m]]$rates[, , s])
    extrapolate <- kannistoExtrapolation(qx, ages.fit, years_chart)
    for (i in 1:nrow(portfolio)){
      ben_s <- ben_s + DFcashflow(extrapolate$qxt, ageStart = portfolio$age[i], omegaAge = 120, 
                                  pensionAge = 65, valyear = valyear, ir = 0.08, type = 1)*pension
      prem_s <- prem_s + DFcashflow(extrapolate$qxt, ageStart = portfolio$age[i], omegaAge = 120, 
                                    pensionAge = 65, valyear = valyear, ir = 0.08, type = 2)*portfolio$Premium[i]
    }
    collectBEL[s, 1] <- ben_s - prem_s
  }
  selectBEL80[[m]] <- quantile(collectBEL, probs = 0.80, type = 1)
  selectBEL90[[m]] <- quantile(collectBEL, probs = 0.90, type = 1)
  selectBEL995[[m]] <- quantile(collectBEL, probs = 0.995, type = 1)
}

SCR80 <- as.numeric(selectBEL80) - BEL
SCR90 <- as.numeric(selectBEL90) - BEL
SCR995 <- as.numeric(selectBEL995) - BEL

df = data.frame(models = c("LC", "RH", "APC", "CBD", "M6", "M7", "M8", "PLAT"), BEL = BEL, SCR80 = SCR80, SCR90 = SCR90, SCR995 = SCR995)
df

# Plot of simulations for Renshaw-Haberman model

pdf("C:/Users/VICTOR/Google Drive/Documents/UoW/MSc Statistics/Dissertation/Dissertation files/graphs/RHsimulate65.pdf", width=7, height=5)
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+horizon))
plot(years,(Dxt/E0xt)["65", ], pch=21, bg='black',
     xlim = range(1991:2055), ylim = range(0.002, 0.04), xlab = "Years", ylab = "Mortality rate at age 65",
     main = "Simulated mortality at age 65")
lines(years_chart, RHqxt["65",], col = "blue", lwd = 2)
matlines(modelSim[[2]]$years, modelSim[[2]]$rates["65", , 1:20], type = "l", lty = 1, col = 3)
abline(h = 0, v = 2014, col = "black")
legend("bottomleft", c("RH", "Simulated paths"), col=c("blue", "green"),
       lty=c(1,1), lwd =2, cex = 0.7, x.intersp=0.3, y.intersp = 0.9, bty = "n")
dev.off()

