library(dplyr)
library(ggplot2)
library(patchwork)

get_NHPP <-function(Tmax, lambda, ...,
                    lambda2=NULL, maxit=10000){
  if(is.null(lambda2)){
    opt <- optimise(lambda, lower = 0, upper = Tmax, maximum = TRUE,
                    ...)
    lambda2 <- opt$objective
  }
  Ts = c()
  Tast = 0
  for(i in 1:maxit){
    Tast = Tast + rexp(1, lambda2)
    if(Tast>Tmax){break}
    if(runif(1) < lambda(Tast, ...)/lambda2){
      Ts =c(Ts,Tast)
    }
  }
  return(Ts)
}


est_intensity <- function(res, Tmax, w){
  ts <- seq(0, Tmax, by=w)
  lamhat <- sapply(1:(length(ts)-1), function(i)sum(ts[i]<res & res<ts[i+1]))/w
  data.frame(time = ts[-1], intensity = lamhat)
}

lambda <- function(x, a, P){a*(cos(x*2*pi/P)+1)}
Lambda <- function(x, a, P){a*((P/(2*pi))*sin(x*2*pi/P)+x)}

rescum <- vector("list", 100)
resint <- vector("list", 100)
set.seed(20260907)
for( i in 1:100 ){
  dat <- get_NHPP(50, lambda, a=2, P=20)
  rescum[[i]] <- data.frame(time=dat, count=seq_along(dat), group=i)
  lamhat <- est_intensity(res = dat, 50, 2.5)
  lamhat$group <- i
  resint[[i]] <- lamhat
}

rescum <- bind_rows(rescum)
resint <- bind_rows(resint)

p1 <- ggplot(rescum, aes(x=time, y=count))+
  geom_step(aes(group = group), linewidth = 0.1)+
  stat_function(fun = Lambda, args = list(a=2, P=20), 
                colour="orangered", linetype=2, linewidth=1)+
  theme_bw(15) + labs(y="cumulative intensity")

p2 <- ggplot(resint, aes(x=time, y=intensity))+
  geom_line(aes(group = group), linewidth = 0.1)+
  stat_function(fun = lambda, args = list(a=2, P=20), 
                colour="orangered", linetype=2, linewidth=1)+
  theme_bw(15)

p1/p2
ggsave("nhpp_cos.png", width = 7, height = 7)

###
loglik <- function(par, y, lambda, Lambda, maxT){
  a <- exp(par[1])
  P <- exp(par[2])
  loglambda <- log(lambda(y, a, P))
  sum(loglambda) - Lambda(maxT, a, P)
}

set.seed(20260908)
dat <- get_NHPP(50, lambda, a=2, P=20)
opt1 <-optim(c(0,3), loglik, y=dat,
             lambda = lambda, Lambda = Lambda,
             maxT = 50,
             control = list(fnscale=-1),
             method = "Nelder-Mead")

ggplot(data = NULL)+
  geom_step(aes(x=dat, y=seq_along(dat))) +
  stat_function(fun = Lambda,
                aes(colour="MLE", linetype="MLE"),
                args = list(a=exp(opt1$par[1]), P=exp(opt1$par[2])), 
                linewidth=1)+
  stat_function(fun = Lambda, 
                aes(colour="true", linetype="true"),
                args = list(a=2, P=20), 
                linewidth=1)+
  scale_linetype_manual(values = 2:3) + 
  labs(colour="", linetype="", x="time", y="cumulative intensity")+
  theme_bw(15)
ggsave("nhpp_fit.png", width = 7, height = 7)
