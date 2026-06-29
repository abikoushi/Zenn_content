library(ggplot2)
library(dplyr)
library(tidyr)
Tmax <- 100
set.seed(1234)
n <- rpois(1,10*Tmax)
ti_c <- sort(runif(n, 0, Tmax))
d <- rbinom(n, 1, 0.2)
ti_d <- sort(ti_c + ifelse(d==1, rlnorm(10000, 2, 0.5), Inf))

Cases <- sapply(1:Tmax, function(t)sum(ti_c<=t))
Deaths <- sapply(1:Tmax, function(t)sum(ti_d<=t))

res_obs <- data.frame(time=1:Tmax, Cases=Cases, Deaths=Deaths) %>% 
  pivot_longer(Cases:Deaths)

ggplot(data = res_obs, aes(x=time, y=value, colour=name, linetype=name))+
  geom_line()+
  theme_bw()
ggsave("obs.png", width = 7, height = 7)
####
res <- vector("list", 100)
Tmax <- 100
set.seed(1111)
for(i in 1:100){
  n <- rpois(1,1000)
  ti_c <- sort(runif(n, 0, Tmax))
  d <- rbinom(n, 1, 0.2)
  ti_d <- sort(ti_c + ifelse(d==1, rlnorm(n, 2, 0.5), Inf))
  
  Cases <- sapply(1:Tmax, function(t)sum(ti_c<=t))
  Deaths <- sapply(1:Tmax, function(t)sum(ti_d<=t))
  
  f <- diff(plnorm(0:Tmax, 2, 0.5))
  
  cases <- diff(c(0,Cases))
  at <- sapply(1:Tmax, function(t){sum(cases[1:t]*f[t:1])})
  At <- cumsum(at)
  
  df <- data.frame(time=1:Tmax,
                   Deaths=Deaths,
                   Cases=Cases,
                   nCFR = Deaths/Cases,
                   cCFR = Deaths/At,
                   group=i)
  
  res[[i]] <- df  
}

res <- bind_rows(res)
res_meth <- pivot_longer(res, nCFR:cCFR, names_to = "method")

ggplot(data = res_meth, aes(x=time, y=value, group=group))+
  geom_line(linewidth=0.1)+
  facet_wrap(~method)+
  geom_hline(yintercept = 0.2, linetype=2)+
  theme_bw()

ggsave("fit.png", width = 7, height = 7)
