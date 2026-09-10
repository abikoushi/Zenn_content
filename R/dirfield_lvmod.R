library(deSolve)
library(ggplot2)
library(tidyr)
library(dplyr)

LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dPrey        <- a*Prey - b*Prey*Predator
    dPredator    <- -c*Predator + d*Prey*Predator
    return(list(c(dPrey, dPredator)))
  })
}

sol_ode_from_states <- function(df_state, times, func, parms){
  res_df <- lapply(1:nrow(df_state), function(i){
    ode_out <- ode(y=c(unlist(df_state[i,,drop=FALSE])), times=times, func=func, parms=parms)
    data.frame(group=i, ode_out)
  })
  dplyr::bind_rows(res_df)
}


pars  <- c(a = 1, 
           b = 2, 
           c = 1,
           d = 3)

yini  <- c(Prey = 1, Predator = 0.5)
times <- seq(0, 50, by = 0.1)
out   <- ode(yini, times, LVmod, pars)


df_out <- pivot_longer(data.frame(out), Prey:Predator) %>% 
  mutate(name=factor(name, levels = c("Prey","Predator")))

ggplot(df_out, aes(x=time, y=value, colour=name, linetype = name))+
  geom_line()+
  theme_bw()

ggsave("lvmod1.png", width = 7, height = 7)

df_state <- expand.grid(Prey=seq(0.1, 2, 0.1),
                        Predator=seq(0.1, 2, 0.1))

times <- seq(0, 0.1, by=0.02)
res <- sol_ode_from_states(df_state=df_state, times=times, func=LVmod, parms=pars)


ggplot(data=res, aes(x=Prey, y=Predator, group=group, color=time)) +
  geom_path(arrow = arrow(length = unit(3, "pt"))) +
  scale_colour_gradient2(low="grey80", high="steelblue")+
  geom_hline(yintercept = pars["a"]/pars["b"], linetype=2)+
  geom_vline(xintercept = pars["c"]/pars["d"], linetype=2)+
  theme_bw()
ggsave("lvmod_dir.png", width = 7, height = 7)
