library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
dat <-read_xls("~/Downloads/y0205000.xls",skip=8) %>% 
  slice(-c(1:8)) %>% 
  select('...1', '人口...3','面積...4') %>% 
  setNames(c("pref","pop","area")) %>% 
  dplyr::filter(!is.na(pop)) %>% 
  mutate(pop=as.integer(pop),area=as.numeric(area)) %>% 
  mutate(pref=gsub("^[0-9][0-9]　", "", pref))

head(dat)

lambahat <- sum(dat$pop)/sum(dat$area)

set.seed(1234)
sim1 <-matrix(rpois(10000*47,dat$area*lambahat),47,10000)
int1 <-as.data.frame(t(apply(sim1,1,quantile,prob=c(0.025,0.975)))) %>% 
  setNames(c("lower","upper")) %>% 
  mutate(area=dat$area)

ggplot(dat,aes(x=area))+
  geom_point(aes(y=pop))+
  geom_ribbon(data=int1,aes(ymin=lower,ymax=upper),alpha=0.3)+
  theme_bw(16)

ggsave("pop_pois.png")

ll <- function(par,y,tau){
  sum(dnbinom(y,size = exp(par[1]), mu = exp(par[2]+tau), log = TRUE))
}

opt <-optim(c(0,0),ll,y=dat$pop,tau=log(dat$area),control = list(fnscale=-1))

r <- rgamma(10000*47,exp(opt$par[1]),exp(opt$par[1]))
sim2 <-matrix(rpois(10000*47,dat$area*r*exp(opt$par[2])),47,10000)
int2 <-as.data.frame(t(apply(sim2,1,quantile,prob=c(0.025,0.975)))) %>% 
  setNames(c("lower","upper")) %>% 
  mutate(area=dat$area)


subdat <- dat[dat$pop > int2$upper,]

ggplot(dat,aes(x=area))+
  geom_point(aes(y=pop))+
  geom_ribbon(data=int2,aes(ymin=lower,ymax=upper),alpha=0.3)+
  geom_text_repel(data = subdat , aes(y=pop, label = pref), family="Osaka")+
  theme_bw(16)

ggsave("pop_negbin.png")
