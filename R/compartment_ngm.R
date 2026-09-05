library(ggplot2)
library(tidyr)
library(dplyr)
library(deSolve)

sol_ode_from_states <- function(df_state, times, func, parms){
  res_df <- lapply(1:nrow(df_state), function(i){
    ode_out <- ode(y=c(unlist(df_state[i,,drop=FALSE])), times=times, func=func, parms=parms)
    data.frame(group=i, ode_out)
  })
  bind_rows(res_df)
}

####
#SIR
####

SIRmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    N <- S + I + R
    dS <- - beta*S*I / N 
    dI <- beta*S*I / N - gamma*I
    dR <- gamma*I 
    return(list(c(dS,  dI, dR)))
  })
}

pars  <- c(beta=0.3, gamma=0.1)
times <- seq(0, 100, by = 0.1)
ini = c(S=999, I=1, R=0)
ode_out <- ode(y=ini, times=times, func=SIRmod, parms=pars)

sir_out <- data.frame(ode_out) |> 
  pivot_longer(S:R) |> 
  mutate(I = name=="I") |> 
  mutate(name = factor(name, levels = c("S","I","R")))

col3 <- hcl.colors(3, palette = "Set2")
ggplot(sir_out, aes(x = time, y = value, colour = name, linetype = name))+
  geom_line() +
  theme_bw() +
  scale_color_manual(values = col3)

ggsave("SIR1.png", width = 6, height = 6)

## exponential growth rate
lambda_sir <- (pars["beta"]-pars["gamma"])
ggplot(sir_out, aes(x=time, y=value))+
  geom_line(aes(group=name,colour=I), show.legend = FALSE)+
  scale_color_manual(values = c("lightgrey",col3[2]))+
  stat_function(fun = function(x)ini[2]*exp(lambda_sir*x),
                linetype=2, colour=col3[2], n=1001, linewidth=1)+
  theme_bw()+labs(y="I") +
  ylim(c(0,1000))

ggsave("SIR_plusexp.png", width = 6, height = 6)

## threshold
df_state <- data.frame(S=seq(1, 991, by=10)) %>% 
  mutate(I=10) %>% 
  mutate(R=1000-(S+I))

res2 <- sol_ode_from_states(df_state=df_state, times=times, func=SIRmod, parms=pars)
res2 <- group_by(res2, group) %>% 
  mutate(iniR = first(S))

ggplot(data=res2, aes(x=time, y=I, colour=iniR, group=group)) +
  geom_line()+
  scale_colour_gradient2(midpoint = 1000*pars["gamma"]/pars["beta"])+
  labs(colour="S(0)")+theme_bw()

ggsave("SIR_thres.png", width = 6, height = 6)

####
#SEIR
####

SEIRmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    N <- S + E + I + R
    dS <- - beta*S*I / N 
    dE <- beta*S*I / N - nu*E
    dI <- nu*E - gamma*I
    dR <- gamma*I 
    return(list(c(dS, dE, dI, dR)))
  })
}

pars  <- c(beta=0.3, nu=0.9, gamma=0.1)
times <- seq(0, 100, by = 0.1)
ini = c(S=999, E=0, I=1, R=0)
ode_out <- ode(y=ini, times=times, func=SEIRmod, parms=pars)

seir_out <- data.frame(ode_out) |> 
  pivot_longer(S:R, names_to = "state") |> 
  mutate(inf = if_else(state=="I"|state=="E", state, "other")) |> 
  mutate(state = factor(state, levels = c("S","E","I","R")))


col4 <- hcl.colors(4, palette = "Set2")

ggplot(seir_out,
       aes(x = time, y = value,
           colour = state, linetype = state))+
  geom_line() +
  theme_bw() +
  scale_color_manual(values = col4)

ggsave("SEIR1.png", width = 6, height = 6)

## exponential growth rate
Solm = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

Sigma <- matrix(c(-pars["nu"], 0,
         pars["nu"], -pars["gamma"]), nrow = 2, byrow = TRUE)

print(solve(Sigma))
print(
  matrix(c(-1/pars["nu"], 0,
           -1/pars["gamma"], -1/pars["gamma"]), nrow = 2, byrow = TRUE)  
)



A <- matrix(c(0, pars["beta"],
              0, 0),nrow = 2, byrow = TRUE) + 
  matrix(c(-pars["nu"], 0,
           pars["nu"], -pars["gamma"]),nrow = 2, byrow = TRUE)
ei_A = eigen(A)
P_A = ei_A$vectors
Pinv_A = solve(ei_A$vectors)

sol1 = sapply(times, Solm, P=P_A, Pinv=Pinv_A, values=ei_A$values, yini=ini[2:3])

ggplot()+
  geom_line(data = seir_out, aes(x=time, y=value, group=state, colour=inf))+
  geom_line(data = NULL,  aes(x=times, y=sol1[1,]), colour=col4[2], linetype = 2, linewidth = 1)+
  geom_line(data = NULL,  aes(x=times, y=sol1[2,]), colour=col4[3], linetype = 2, linewidth = 1)+
  theme_bw() + 
  scale_colour_manual(values = c(col4[2:3], "grey"))+
  ylim(c(0,1050))+labs(colour = "state")

ggsave("SEIR_plusexp.png", width = 6, height = 6)

## threshold
df_state <- data.frame(S=seq(1, 991, by=10)) %>% 
  mutate(E=0,I=10) %>% 
  mutate(R=1000-(S+E+I))

res2 <- sol_ode_from_states(df_state=df_state, times=times, func=SEIRmod, parms=pars)
res2 <- group_by(res2, group) %>% 
  mutate(iniR = first(S))

ggplot(data=res2, aes(x=time, y=E+I, colour=iniR, group=group)) +
  geom_line()+
  scale_colour_gradient2(midpoint = 1000*pars["gamma"]/pars["beta"])+
  labs(colour="S(0)")+theme_bw()

ggsave("SEIR_thres.png", width = 6, height = 6)
