library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(dplyr)
###
mod <- cmdstan_model("./stan/instrument_stan/IV.stan")
rand_pre = function(n, alpha_xu, alpha_xz, alpha_yx, alpha_yu){
  U <- rnorm(n)
  Z <- rnorm(n)
  X <- rnorm(n, alpha_xz*Z+alpha_xu*U)
  Y <- rnorm(n, alpha_yx*X+alpha_yu*U)
  return(data.frame(Y,Z,X))
}

rand_post_average = function(n,alpha_xu, alpha_xz, alpha_yx, alpha_yu){
  U <- rnorm(n)
  
  ##X=1
  Y1 <- rnorm(n, alpha_yx+alpha_yu*U)
  
  ##X=0
  Y0 <- rnorm(n, alpha_yu*U)
  return(mean(Y1-Y0))
}

alpha_xu = 0.5
alpha_xz = 0.9
alpha_yx = 0.5
alpha_yu = 0.2

#ATE_true <- rand_post_average(100000, alpha_xu, alpha_xz, alpha_yx, alpha_yu) 
#print(ATE_true) #0.5

set.seed(1234)
dat <- rand_pre(500, alpha_xu, alpha_xz, alpha_yx, alpha_yu)
datlist <- as.list(dat)
datlist$N <- nrow(dat)

fit <- mod$sample(
  data = datlist,
  seed = 12345,
  chains = 4,
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

mcmc_trace(fit$draws(c("alpha_xu", "alpha_xz", "alpha_yx", "alpha_yu")))

ggsave("trace.png", width = 7, height = 7)

TE <- reshape2::melt(fit$draws("TE"))
ATE <- summarise(TE, ATE = mean(value)) %>% 
  pull(ATE)

p1 <- ggplot(TE)+
  geom_histogram(aes(x=value), bins=40)+
  geom_vline(xintercept = alpha_yx, colour="grey", linetype=2)+
  annotate(geom = "point", x=ATE, y=-1 , colour="steelblue", size=4, shape=17) +
  theme_classic()+labs(x="treatment effect")

print(p1)
ggsave("hist1.png", plot = p1, width = 5, height = 5)

hatalpha_yx = with(dat, as.vector(solve(t(Z)%*%X, t(Z)%*%Y)))
hatalpha_yx

p2 <- p1 + geom_vline(xintercept = alpha_yx, colour="orange", linetype=3)

ggsave("hist2.png", plot = p2, width = 5, height = 5)
