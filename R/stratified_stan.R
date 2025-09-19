#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)
library(posterior)
library(bayesplot)
library(patchwork)
library(cmdstanr)

df_prob = expand.grid(y=0:1, x=0:1, z=0:1) %>%
  mutate(p = c(0.2,0.05,0.2,0.2,0.05,0.245,0.005,0.05))
kable(df_prob, digits = 3)


rand_case1 = function(n, xi, psi, gamma){
  ## draw x given by z
  z = rbinom(n, 1, gamma)
  x = rbinom(n, 1, psi[z+1])
  y = rbinom(n, 1, xi[cbind(z+1,x+1)])
  return(data.frame(Y=y, X=x, Z=z))
}

rand_case2 = function(n, xi, phi, delta){
  ## draw z given by x
  x = rbinom(n, 1, delta)
  z = rbinom(n, 1, phi[x+1])
  y = rbinom(n, 1, xi[cbind(z+1,x+1)])  
  return(data.frame(Y=y,X=x,Z=z))
}

rand_ATE_case1 = function(n, xi, psi, gamma){
  z = rbinom(n, 1, gamma)
  y0 = rbinom(n, 1, xi[cbind(z+1,1)])
  y1 = rbinom(n, 1, xi[cbind(z+1,2)])
  return(y1-y0)
}

rand_ATE_case2 = function(n, xi, phi, delta){
  z0 = rbinom(n, 1, phi[1])
  y0 = rbinom(n, 1, xi[cbind(z0+1,1)])
  z1 = rbinom(n, 1, phi[2])
  y1 = rbinom(n, 1, xi[cbind(z1+1,2)])  
  return(y1-y0)
}

ATE_case1 = function(xi, psi, gamma){
  y0 = (1-gamma)*xi[1,1]+gamma*xi[2,1]
  y1 = (1-gamma)*xi[1,2]+gamma*xi[2,2]
  return(y1-y0)
}

rand_ATE_case2 = function(n, xi, phi, delta){
  z0 = rbinom(n, 1, phi[1])
  y0 = rbinom(n, 1, xi[cbind(z0+1,1)])
  z1 = rbinom(n, 1, phi[2])
  y1 = rbinom(n, 1, xi[cbind(z1+1,2)])  
  return(y1-y0)
}

ATE_case2 = function(xi,  phi, delta){
  y0 = phi[1]*xi[2,1] + (1-phi[1])*xi[1,1]
  y1 = phi[2]*xi[2,2] + (1-phi[2])*xi[1,2]
  return(y1-y0)
}

loglik_case1 = function(y, z, x, xi, psi, gamma){
  # x given by z
  sum(dbinom(z, 1, gamma, log = TRUE))+
    sum(dbinom(x, 1, psi[z+1], log = TRUE))+
    sum(dbinom(y, 1, xi[cbind(z+1,x+1)], log = TRUE))  
}

loglik_case2 = function(y, z, x, xi, phi, delta){
# z given by x
sum(dbinom(x, 1, delta, log = TRUE)) + 
  sum(dbinom(z, 1, phi[x+1], log = TRUE)) + 
  sum(dbinom(y, 1, xi[cbind(z+1,x+1)], log = TRUE))
}

PY_x = group_by(df_prob, x) %>% 
  summarise(p = sum(y*p)/sum(p))

PY_xz = group_by(df_prob, x, z) %>% 
  summarise(p = sum(y*p)/sum(p))

PX_z = group_by(df_prob,  z) %>% 
  summarise(p = sum(x*p)/sum(p))

PZ_x = group_by(df_prob,  x) %>% 
  summarise(p = sum(z*p)/sum(p))

PX = summarise(df_prob, p = sum(x*p)/sum(p))
PZ = summarise(df_prob, p = sum(z*p)/sum(p))

Xi = matrix(PY_xz$p, 2, 2)
rownames(Xi) = paste0("z", 0:1)
colnames(Xi) = paste0("x", 0:1)

dat1 = rand_case1(n = 10000, xi = Xi, psi = PX_z$p, gamma = PZ$p)
tab1 = group_by(dat1,Y,X,Z) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(p=n/sum(n)) %>% 
  arrange(Z,X,Y)


dat2 = rand_case2(n = 10000, xi = Xi, phi = PZ_x$p, delta = PX$p)
tab2 = group_by(dat2,Y,X,Z) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(p=n/sum(n)) %>% 
  arrange(Z,X,Y)

kable(cbind(df_prob, "rand_case1" = tab1$p, "rand2_case"=tab2$p),
      digits=3)

loglik_case1(y = dat1$Y, z = dat1$Z, x=dat1$Z,
             xi = Xi, psi = PX_z$p, gamma = PZ$p)
loglik_case2(y = dat1$Y, z = dat1$Z, x=dat1$Z,
             xi = Xi, phi = PZ_x$p, delta = PX$p)

loglik_case1(y = dat2$Y, z = dat2$Z, x=dat2$Z,
             xi = Xi, psi = PX_z$p, gamma = PZ$p)
loglik_case2(y = dat2$Y, z = dat2$Z, x=dat2$Z,
             xi = Xi, phi = PZ_x$p, delta = PX$p)


###
check_cmdstan_toolchain()
PATH_TO_CMDSTAN = cmdstan_path()
set_cmdstan_path(PATH_TO_CMDSTAN)

mod1 <- cmdstan_model("./stan/stratified_stan/case1.stan")
mod2 <- cmdstan_model("./stan/stratified_stan/case2.stan")

set.seed(1234)
dat = rand_case1(n = 500, xi = Xi, psi = PX_z$p, gamma = PZ$p)
# names correspond to the data block in the Stan program
data_list <- list(N = nrow(dat),
                  Y = dat$Y, X=dat$X, Z=dat$Z,
                  alpha=1)

fit1 <- mod1$sample(
  data = data_list,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2500,
  refresh = 500 # print update every 500 iters
)

fit2 <- mod2$sample(
  data = data_list,
  seed = 2,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2500,
  refresh = 500 # print update every 500 iters
)


parnames = c("lp__","Xi","psi","phi","delta","gamma")
p_tr1 = mcmc_trace(fit1$draws(parnames))+
  ggtitle("case 1")+theme_bw()
p_tr2 = mcmc_trace(fit2$draws(parnames))+
  ggtitle("case 2")+theme_bw()
p_tr = p_tr1 / p_tr2
print(p_tr)
ggsave(filename = "traceplot.png", plot = p_tr, width = 10, height = 8)

df_par = bind_rows(
  mutate(reshape2::melt(fit1$draws(parnames)), model=1), 
  mutate(reshape2::melt(fit2$draws(parnames)), model=2) ) %>% 
  mutate(model = factor(model))

p_ecdf = ggplot(df_par, aes(x=value, colour = model, linetype = model))+
  stat_ecdf()+
  facet_wrap(~variable, scales = "free_x", nrow=3)+
  theme_bw(16)
print(p_ecdf)
ggsave(filename = "param_ecdf.png", plot = p_ecdf, width = 14, height = 8)

fit1$summary()
fit2$summary()

#Xihat = matrix(fit1$summary("Xi", c("mean","sd"))$mean, 2, 2)

summary1 = fit1$summary(c("Xi","psi","phi","delta","gamma"), c("mean","sd"))
summary2 = fit2$summary(c("Xi","psi","phi","delta","gamma"), c("mean","sd"))

cbind(summary1[,1], model1=summary1[,-1], model2=summary2[,-1],
      true=c(PY_xz$p,PX_z$p,PZ_x$p,PZ$p,PX$p)) %>% kable(digits = 2)

rbind(
  mutate(fit1$summary(c("D1","D2"),c("mean","sd")), variable=paste(variable,"(model1)")),
  mutate(fit2$summary(c("D1","D2"),c("mean","sd")), variable=paste(variable,"(model2)"))  
) %>% kable(digits = 2)
  

print(round(ATE_case1(xi = Xi, psi=PX_z$p, gamma=PZ$p),3))
print(round(ATE_case2(xi = Xi, phi=PZ_x$p, delta=PX$p),3))
