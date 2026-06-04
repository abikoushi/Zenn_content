library(ggplot2)
library(deSolve)
library(dplyr)
library(tidyr)

SIRmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS <- - beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I 
    return(list(c(dS, dI, dR)))
  })
}

pars_sir  <- c(beta=1, gamma=0.25)
ini_sir = c(S=0.999, I=0.001, R=0)
times <- seq(0, 50, by = 0.1)
ode_out <- ode(y=ini_sir, times=times, func=SIRmod, parms=pars_sir)

sir_out <- data.frame(ode_out) %>% 
  pivot_longer(S:R) %>% 
  mutate(I = name=="I")

lambda_sir <- (pars_sir["beta"]-pars_sir["gamma"])

p1 <- ggplot(sir_out, aes(x=time, y=value))+
  geom_line(aes(group=name,colour=I), show.legend = FALSE)+
  scale_color_manual(values = c("lightgrey","black"))+
  stat_function(fun = function(x)ini_sir[2]*exp(lambda_sir*x),
                linetype=2, colour="royalblue", n=1001, linewidth=1)+
  theme_bw()+labs(y="I") +
  ylim(c(0,1.05))

print(p1)

ggsave("SIR_expgrowth.png", width = 7, height = 7)

####

R0hat_gamma <- function(m,CV){
  r <- 0.135
  ab2 <- (m*CV)^2
  b <- ab2/m
  a <- m/b
  int <- integrate(function(x){exp(-r*x)*dgamma(x,shape = a,scale = b)},0,Inf)
  return(1/int$value)  
}

resdf <- expand.grid(mean=seq(1,6,by=0.1),
                     CV=seq(0.1,1,by=0.1)) %>% 
  rowwise() %>% 
  mutate(R0=R0hat_gamma(mean, CV))


ggplot(resdf, aes(x=mean, y=CV, fill=R0))+
  geom_tile(width = 0.1, height = 0.1)+
  scale_fill_viridis_c()

ggsave("fig3_3.png", width = 7, height = 7)
