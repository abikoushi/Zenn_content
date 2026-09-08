library(ggplot2)
library(tidyr)
library(dplyr)
library(deSolve)
library(nleqslv)

SIRmod <- function(Time, State, Pars) {
  n <- Pars$n
  beta <- Pars$beta
  gamma <- Pars$gamma
  N <- Pars$N
  S <- State[1:n]
  I <- State[(n+1):(2*n)]
  R <- State[(2*n+1):(3*n)]
  dS <- - drop(beta%*%I)*S/N
  dI <- drop(beta%*%I)*S/N - gamma*I
  dR <- gamma*I
  return(list(c(dS, dI, dR)))
}

set.seed(20260909); beta <- matrix(rexp(4,1), byrow = TRUE, nrow = 2)
print(beta)

N <- c(1, 1)
pars  <- list(beta = beta, gamma = 0.1, N=N, n=2)
times <- seq(0, 100, by = 0.1)

ini = c(S1=1, S2=0.999,
        I1=0, I2=0.001,
        R1=0, R2=0)

ode_out <- ode(y=ini, times=times, func=SIRmod, parms=pars)

sir_out <- data.frame(ode_out) |> 
  pivot_longer(S1:R2) |> 
  mutate(state = substr(name, 1,1), 
         subgroup = substr(name, 2,2)) |> 
  mutate(state = factor(state, levels = c("S","I","R")))

col2 <- hcl.colors(2, palette = "Set2")
ggplot(sir_out, aes(x = time, y = value,
                    colour = subgroup, 
                    linetype = subgroup))+
  geom_line() + 
  facet_grid(row=vars(state)) + 
  scale_color_manual(values = col2)+
  theme_bw(15)

ggsave("SIR1.png", width = 7, height = 7)

####
#exp growth
####
A <- pars$beta + diag(-pars$gamma, pars$k)
eiA <- eigen(A)

matexp = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

ei_A = eigen(A)
P_A = ei_A$vectors
Pinv_A = solve(ei_A$vectors)

sol0 = sapply(times, matexp, 
              P=P_A, Pinv=Pinv_A, values=ei_A$values,
              yini=ini[3:4])
row.names(sol0) <- names(ini[3:4])
exp_out <- data.frame(time=times, t(sol0)) |> 
  pivot_longer(I1:I2)
sir_out_inf <- dplyr::filter(sir_out, state=="I")

ggplot(sir_out_inf, mapping = aes(x=time, y=value, group=name, colour=name))+
  geom_line()+
  geom_line(data = exp_out, linetype = 2, linewidth = 1)+
  scale_colour_manual(values = col2)+
  labs(colour = "state") +
  facet_grid(row=vars(name))+
  ylim(c(0,1.05)) +
  theme_bw(15)

ggsave("SIR_exp.png", width = 7, height = 7)

#####
#final size
#####

KL <- pars$beta%*%diag(1/pars$gamma, pars$k)
ei_KL <- eigen(KL)

sir_out_rem <- dplyr::filter(sir_out, grepl("^R",name))

fz <- function(z){z + expm1(-KL%*%z)}
res <- nleqslv(c(1,1), fz)
df_z <- data.frame(z = res$x, 
                   name = paste0("R",1:2))

ggplot(sir_out_rem, aes(x = time, y = value, colour = subgroup))+
  geom_hline(data = df_z,  aes(yintercept = z), linetype = 2)+
  geom_line() + 
  facet_grid(row=vars(name)) + 
  scale_color_manual(values = col2)+
  theme_bw(15)+ylim(c(0,1))

ggsave("SIR_R.png", width = 7, height = 7)
