
DFcashflow = function(qxt, ageStart, omegaAge, pensionAge, valyear, ir, type = 1){
  if (type == 1){
    CF <- c(rep(0,max(0,(pensionAge-ageStart+1))), rep(1,(omegaAge-max(pensionAge, ageStart-1))))
  }else{
    CF <- c(rep(1,max(0,(pensionAge-ageStart+1))), rep(0,(omegaAge-max(pensionAge, ageStart-1))))
  }
  
  t_start <- which(colnames(qxt) == toString(valyear))
  x_start <- which(rownames(qxt) == toString(ageStart))
  qx      <- diag(qxt[x_start:(x_start+(omegaAge-ageStart)), t_start:(t_start+(omegaAge-ageStart))])
  pxt <- cumprod(1- qx)
  k <- 0:(omegaAge-ageStart)
  interestRates <- rep(ir, (omegaAge-ageStart)+1)
  v <- (1 + interestRates)^-k
  dfCF <- sum(pxt*v*CF)
  dfCF
}

pxt <- function(qxt, valyear, ageStart, omegaAge){
  t_start <- which(colnames(qxt) == toString(valyear))
  x_start <- which(rownames(qxt) == toString(ageStart))
  qxt <- diag(qxt[x_start:(x_start+(omegaAge-ageStart)), t_start:(t_start+(omegaAge-ageStart))])
  pxt <- c(1, head(cumprod(1- qxt), -1))
  pxt
}

kannistoExtrapolation = function(qx, ages, years, max_age=120, nObs=15){
  agesExtrap <- (tail(ages,1)+1):max_age
  
  obsAges <- tail(ages, nObs)

  extrapolateQxt <- sapply(years, function(x){
    obsQx <- pmin(tail(qx[,toString(x)], nObs),1)
    obsQxt  <- obsQx 
    obsQx <- -log(1 - obsQx)
    obsQx <- log(obsQx/(1 - obsQx))
    
    if (length(na.omit(obsQx))>2){
      model <- lm(formula = obsQx~obsAges)
      phi1 <- exp(model$coefficient[1])
      phi2 <- model$coefficient[2]
      1-exp(-phi1*exp(phi2*agesExtrap)/(1+phi1*exp(phi2*agesExtrap)))
    }else{
      seq(obsQxt[nObs], to = 1, length = 30)
    }
  })
  
  qxt <- rbind(qx,  extrapolateQxt)
  ages <- c(ages, agesExtrap)
  
  list(qxt = qxt, ages.fit = ages)
}