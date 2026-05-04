library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

SIRmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS <- - beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I 
    return(list(c(dS,  dI, dR)))
  })
}

sol_ode_from_states <- function(df_state, times, func, parms){
  res_df <- lapply(1:nrow(df_state), function(i){
    ode_out <- ode(y=c(unlist(df_state[i,,drop=FALSE])), times=times, func=func, parms=parms)
    data.frame(group=i, ode_out)
  })
  bind_rows(res_df)
}


pars  <- c(beta=1.5, gamma=1)
times <- seq(0, 20, by = 0.1)
ini = c(S=0.99, I=0.01, R=0)
ode_out1 <- ode(y=ini, times=times, func=SIRmod, parms=pars)

R0 <- pars["beta"]/pars["gamma"]


z_root <- uniroot(function(z)R0 + (log1p(-z)-log(ini["S"]))/(z-ini["R"]), interval = c(0.001,1))

z_root$root

res1 <- data.frame(ode_out1) %>% 
  pivot_longer(S:R) %>% 
  mutate(name=factor(name, levels = c("S","I","R")))

ggplot(data=res1, aes(x=time, y=value, colour=name, linetype=name)) +
  geom_hline(yintercept = z_root$root, linetype=2)+
  geom_line()+
  theme_bw()

ggsave("SIR_sol.png", width = 5, height = 5)

####

df_state <- data.frame(S=seq(0.01, 0.99, 0.01)) %>% 
  mutate(I=0.01) %>% 
  mutate(R=1-(S+I))

head(df_state)

times <- seq(0, 20, by = 0.1)
res2 <- sol_ode_from_states(df_state=df_state, times=times, func=SIRmod, parms=pars)
res2 <- group_by(res2, group) %>% 
  mutate(iniR = first(S))

ggplot(data=res2, aes(x=time, y=I, colour=iniR, group=group)) +
  geom_line()+
  #scale_color_viridis_c()+
  scale_colour_gradient2(midpoint = pars["gamma"]/pars["beta"])+
  labs(colour="S(0)")+theme_bw()
ggsave("SIR_I.png", width = 5, height = 5)

###

zfun <- function(R0){
  z_root <- uniroot(function(z){R0 + log1p(-z)/z}, interval = c(1e-8,1))
  z_root$root
}
rv <-seq(1.1,10,by=0.1)
zv <- sapply(rv, zfun)
ggplot(data = NULL)+
  geom_line(aes(x=rv, y=zv)) + 
  theme_bw()+labs(x="R0", y="z")

ggsave("size_by_R0.png", width = 5, height = 5)
