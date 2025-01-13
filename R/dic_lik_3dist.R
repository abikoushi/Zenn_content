#https://github.com/abikoushi/IntervalCensored.jl/blob/main/src/calclp.jl
#install.packages("coarseDataTools")
library(coarseDataTools)
library(dplyr)
library(ggplot2)
library(rbenchmark)
#gamma dist
eqgamma <- function(x, shape, scale){
  pgamma(x/scale, shape+1)+(x/scale)*pgamma(x/scale, shape, lower.tail=FALSE)/shape
}
meangamma <- function(shape, scale){
  scale*shape
}

#weibull dist
eqweibull <- function(x, shape, scale){
  pgamma((x/scale)^(shape), 1/shape)
}
meanweibull <- function(shape, scale){
  scale*gamma(1+1/shape)
}


#lognormal dist
erf <- function(x){
  2*pnorm(sqrt(2)*x) - 1
}

erfc <- function(x){
  2*pnorm(-sqrt(2)*x) 
}
eqlnorm <- function(x, meanlog, sdlog){
  ifelse(
    x < 0, 
    0,
    0.5*(1+x*erfc((log(x) - meanlog)/(sdlog*sqrt(2)))*exp(-(meanlog + (sdlog^2)/2)) - erf((meanlog + sdlog^2 - log(x))/(sdlog * sqrt(2))))
  )
}
meanlnorm <- function(mu, sigma){
  exp(mu+sigma^2/2)
}

##
#check 3 dist
#weibull
eqweibull0 <- function(y, shp, scl){
  integrate(function(x)pweibull(x, shp, scale = scl, lower.tail = FALSE)/meanweibull(shp, scl), 0, y)$value
}
eqweibull0 <- Vectorize(eqweibull0)
curve(eqweibull0(x, 1.3, 1.5), 0, 10)
curve(eqweibull(x, 1.3, 1.5), add=TRUE, col="orangered", lty=2)

#gamma
eqgamma0 <- function(y, shp, scl){
  integrate(function(x)pgamma(x, shp, scale = scl, lower.tail = FALSE)/meangamma(shp, scl), 0, y)$value
}
eqgamma0 <- Vectorize(eqgamma0)
curve(eqgamma0(x, 1.9, 1.5), 0, 10)
curve(eqgamma(x, 1.9, 1.5), add=TRUE, col="orangered", lty=2)

#lognormal
eqlnorm0 <- function(y, mu, sigma){
  integrate(function(x)plnorm(x, mu, sigma, lower.tail = FALSE)/meanlnorm(mu, sigma), 0, y)$value
}
eqlnorm0 <- Vectorize(eqlnorm0)
curve(eqlnorm0(x, 1.1, 0.5), 0, 10)
curve(eqlnorm(x, 1.1, 0.5), add=TRUE, col="orangered", lty=2)
##

#doubly interval censored
evaluatelp0 <- function(EL, ER, SL, SR, par, eqcdf, meanf){
  ll = 0 
  mu = meanf(par)
  if(ER <= SL){
    ll = log(eqcdf(SR-ER, par) - eqcdf(SL-ER, par) - (eqcdf(SR-EL, par) - eqcdf(SL-EL, par))) + log(mu)
  }else if (SL < ER & ER < SR){
    ll = log((ER - SL) + mu*(eqcdf(SR-ER, par) - (eqcdf(SR-EL, par) - eqcdf(SL-EL, par))))
  }else if(SR <= ER){
    ll = log((SR - SL) - mu*(eqcdf(SR-EL, par) - eqcdf(SL-EL, par)) )  
  }
  return(ll)
}

# single interval
evaluatelp1 <- function(EL, ER, SL, SR, par, cdf){
  ll = log(cdf(SR-EL, par) - cdf(SL-EL, par))
  return(ll)
}

# no censored
evaluatelp2 <- function(EL, ER, SL, SR, par, logpdf){
  logpdf(SL-EL, par)
}

eval_lp <- function(par, EL, ER, SL, SR, type, dist){
  if(dist=="weibull"){
    eqcdf <- function(x, par){
      eqweibull(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    cdf <- function(x, par){
      pweibull(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    logpdf<- function(x, par){
      dweibull(x, shape = exp(par[1]), scale = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meanweibull(exp(par[1]), exp(par[2]))
    }
  }else if(dist == "gamma"){
    eqcdf <- function(x, par){
      eqgamma(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    cdf <- function(x, par){
      pgamma(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    logpdf<- function(x, par){
      dgamma(x, shape = exp(par[1]), scale = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meangamma(exp(par[1]), exp(par[2]))
    }
  }else if(dist == "lognormal"){
    eqcdf <- function(x, par){
      eqlnorm(x, meanlog = par[1], sdlog = exp(par[2]))
    }
    cdf <- function(x, par){
      plnorm(x,  meanlog = par[1], sdlog = exp(par[2]))
    }
    logpdf<- function(x, par){
      dlnorm(x, meanlog = par[1], sdlog = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meanlnorm(par[1], exp(par[2]))
    }
  }else{
    warning("this distribution is not implemented")
  }
  N <- length(EL)
  ll = 0
  for(i in 1:N){
    if(type[i]==0){
      ll = ll + evaluatelp0(EL[i], ER[i], SL[i], SR[i], par, eqcdf = eqcdf, meanf = meanf)
    }else if(type[i]==1){
      ll = ll + evaluatelp1(EL[i], ER[i], SL[i], SR[i], par, cdf = cdf)
    }else if(type[i]==2){
      ll = ll + evaluatelp2(EL[i], ER[i], SL[i], SR[i], par, logpdf = logpdf)
    }else{
      warning("this censor type is not implemented")
    }    
  }
  return(ll)
}

data(fluA.inc.per)
head(fluA.inc.per)


set.seed(12345); inip <- rnorm(2)
coarseDataTools::loglikhd(inip, dat=fluA.inc.per, dist = "W")
with(fluA.inc.per, eval_lp(inip, EL, ER, SL, SR, type = type, dist = "weibull"))

coarseDataTools::loglikhd(inip, dat=fluA.inc.per, dist = "G")
with(fluA.inc.per, eval_lp(inip, EL, ER, SL, SR, type = type, dist = "gamma"))

coarseDataTools::loglikhd(inip, dat=fluA.inc.per, dist = "L")
with(fluA.inc.per, eval_lp(inip, EL, ER, SL, SR, type = type, dist = "lognormal"))


lpW <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "weibull"))  
}

lpG <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "gamma")) 
}

lpL <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "lognormal")) 
}

bm_w <- benchmark(flu_opt_w0 <- optim(c(0,1), loglikhd, dat=fluA.inc.per, dist = "W"),
                  flu_opt_w <- optim(c(0,1), lpW, dat=fluA.inc.per), order = NULL)



bm_g <- benchmark(flu_opt_g0 <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "G"),
                  flu_opt_g <- optim(c(0,1), lpG, dat=fluA.inc.per), order = NULL)


bm_l <- benchmark(flu_opt_l0 <- optim(c(0,1),loglikhd, dat=fluA.inc.per,dist = "L"),
                  flu_opt_l <- optim(c(0,1), lpL, dat=fluA.inc.per), order = NULL)


print(bm_g)
print(all.equal(flu_opt_g0$par, flu_opt_g$par))

print(bm_w)
print(all.equal(flu_opt_w0$par, flu_opt_w$par))

print(bm_l)
print(all.equal(flu_opt_l0$par, flu_opt_l$par))

#png("fluA.png")
c3 <- c("royalblue", "orangered", "forestgreen")
curve(pweibull(x, exp(flu_opt_w$par[1]), exp(flu_opt_w$par[2]), lower.tail = FALSE),
      xlim=c(0,5), col=c3[1], lty=2, ylab="CCDF", xlab="incubation period")
curve(pgamma(x, exp(flu_opt_g$par[1]), scale = exp(flu_opt_g$par[2]), lower.tail = FALSE),
      add=TRUE, col=c3[2], lty=3)
curve(plnorm(x, flu_opt_l$par[1], sdlog = exp(flu_opt_l$par[2]), lower.tail = FALSE),
      add=TRUE, col=c3[3], lty=4)
legend("topright",  legend = c("Weibull", "gamma", "lognormal"), lty=2:4, col=c3, lwd=1)
#dev.off()
