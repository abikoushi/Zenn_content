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

SIRmod <- function(Time, State, Pars) {
  k <- Pars$k
  beta <- Pars$beta
  gamma <- Pars$gamma
  S <- State[1:k]
  I <- State[(k+1):(2*k)]
  R <- State[(2*k+1):(3*k)]
  dS <- - drop(beta%*%I)*S
  dI <- drop(beta%*%I)*S - gamma*I
  dR <- gamma*I
  return(list(c(dS, dI, dR)))
}

# beta <- matrix(c(0.3, 0.2, 0.2,
#                  0.2, 0.2, 0.1,
#                  0.2, 0.1, 0.1), byrow = TRUE, nrow = 3)

set.seed(20260908); beta <- matrix(rexp(9,10), byrow = TRUE, nrow = 3)
print(beta)

pars  <- list(beta = beta, gamma = 0.1, k=3)
times <- seq(0, 150, by = 0.1)

ini = c(S1=1, S2=0.999, S3=1,
        I1=0, I2=0.001, I3=0,
        R1=0, R2=0, R3=0)

ode_out <- ode(y=ini, times=times, func=SIRmod, parms=pars)


sir_out <- data.frame(ode_out) |> 
  pivot_longer(S1:R3) |> 
  mutate(state = substr(name, 1,1), 
         subgroup = substr(name, 2,2)) |> 
  mutate(state = factor(state, levels = c("S","I","R")))

col3 <- hcl.colors(3, palette = "Set2")
ggplot(sir_out, aes(x = time, y = value,
                    colour = subgroup, 
                    linetype = subgroup))+
  geom_line() + 
  facet_grid(row=vars(state)) + 
  scale_color_manual(values = col3)+
  theme_bw(15)

####
#exp growth
####
A <- pars$beta + diag(-pars$gamma,3)
eiA <- eigen(A)

matexp = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

ei_A = eigen(A)
P_A = ei_A$vectors
Pinv_A = solve(ei_A$vectors)

sol0 = sapply(times, matexp, 
              P=P_A, Pinv=Pinv_A, values=ei_A$values,
              yini=ini[4:6])
row.names(sol0) <- names(ini[4:6])
exp_out <- data.frame(time=times, t(sol0)) |> 
  pivot_longer(I1:I3)
sir_out_inf <- dplyr::filter(sir_out, state=="I")

ggplot(sir_out_inf, mapping = aes(x=time, y=value, group=name, colour=name))+
  geom_line()+
  geom_line(data = exp_out, linetype = 2, linewidth = 1)+
  scale_colour_manual(values = c(col3, "grey"))+
  labs(colour = "state") +
  facet_grid(row=vars(name))+
  ylim(c(0,1.05)) +
  theme_bw(15)


####
#threshold
####

R <- pars$beta%*%diag(1/pars$gamma,3)
ei_R <- eigen(R)
ei_R$values

df_state <- data.frame(S1=seq(0.01, 0.99, by=0.02)) |> 
  mutate(I1=0.01) |> 
  mutate(R1=1-(S1+I1)) |>
  mutate(S2=S1,S3=S1,
         I2=0,I3=0) |> 
  mutate(R2=1-(S2+I2),
         R3=1-(S3+I3)) |>
  dplyr::select(S1,S2,S3,I1,I2,I3,R1,R2,R3)


res1 <- sol_ode_from_states(df_state=df_state, 
                            times = times,
                            func=SIRmod, parms=pars)

res1 <- group_by(res1, group) |>  
  mutate(iniS1 = first(S1))

ggplot(data=res1, aes(x=time, y=I1+I2+I3, group = group, colour=iniS1)) +
  geom_line()+
  scale_colour_viridis_c()+
  theme_bw() + labs(colour="S(0)")


res2 <- sol_ode_from_states(df_state=df_state, 
                            times = seq(0, 1, by=0.01),
                            func=SIRmod, parms=pars)

res2 <- group_by(res2, group) |>  
  mutate(iniS1 = first(S1))

ggplot(data=res2, aes(x=iniS1, y=I1+I2+I3, group = group, colour=time)) +
  geom_line(arrow = arrow(length = unit(3,"pt")))+
  geom_vline(xintercept =  1/ei_R$values[1], linetype = 2)+
  scale_colour_gradient2(low="grey80", high="steelblue")+
  labs(x="S(0)", y="I(1)")+
  theme_bw()


#####
#final size
#####

sir_out_rem <- dplyr::filter(sir_out, grepl("^R",name))

ggplot(sir_out_rem, aes(x = time, y = value, colour = subgroup))+
  geom_line() + 
  facet_grid(row=vars(name)) + 
  scale_color_manual(values = col3)+
  theme_bw(15)+ylim(c(0,1))


optim(c(0,0,0), function(z){sum(abs(-expm1(-R%*%z)))})
