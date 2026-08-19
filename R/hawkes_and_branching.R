library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

simulate_branching_death <- function(mu,
                                     alpha,
                                     beta,
                                     Tmax) {
  
  N0 <- rpois(1, Tmax * mu)
  if(N0==0L){
    return(
      data.frame(
        id = NA,
        parent = 0L,
        birth_time = NA,
        death_time = NA,
        generation = NA
      )
    )
  }
  
  t0 <- sort(runif(N0,0,Tmax))
  individuals <- data.frame(
    id = seq_len(N0),
    parent = 0L,
    birth_time = t0,
    death_time = t0 + rexp(N0, rate = beta),
    generation = 1L
  )
  
  next_id <- N0 + 1L
  i <- 1L
  
  while (i <= nrow(individuals)) {
    
    birth <- individuals$birth_time[i]
    death <- min(individuals$death_time[i], Tmax)
    
    if (birth < death) {
      
      t <- birth
      
      repeat {
        
        t <- t + rexp(1, rate = alpha)
        
        if (t > death || t > Tmax) {
          break
        }
        
        individuals <- rbind(
          individuals,
          data.frame(
            id = next_id,
            parent = individuals$id[i],
            birth_time = t,
            death_time = t + rexp(1, rate = beta),
            generation =
              individuals$generation[i] + 1L
          )
        )
        
        next_id <- next_id + 1L
      }
    }
    
    i <- i + 1L
  }
  return(individuals[order(individuals$birth_time),])
}


Tau <- 10
mu <- 0.5
a <- 0.8
b <- 1
set.seed(20260819)
dat <- simulate_branching_death(mu, a, b, Tau)

segments <- dplyr::filter(dat, parent > 0) |> 
  rowwise() |> 
  mutate(
    x = dat$birth_time[dat$id == parent],
    y = dat$generation[dat$id == parent],
    xend = birth_time,
    yend = generation-0.05
  )

ggplot(dat) +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.02, "npc")),
    colour="grey") +
  geom_segment(aes(x = birth_time, xend=death_time, y = generation, yend=generation),
               position = position_jitter(height = 0.05, seed = 123))+
  geom_point(aes(birth_time, generation), size=3, shape=1,
             position = position_jitter(height = 0.05, seed = 123)) +
  geom_point(aes(death_time, generation), size=3, shape=4,
             position = position_jitter(height = 0.05, seed = 123)) +
  geom_rug(aes(x=birth_time))+
  scale_y_reverse() +
  labs(x = "time", y="generation") +
  theme_bw()

ggsave("birth_death.png", width = 7, height = 7)

xv <- seq(0,10,by=0.01)

bdmf <- data.frame(
  time = xv,
  Birth = sapply(xv, function(t){
    sum(dat$birth_time <t)
  }),
  Current = sapply(xv, function(t){
    sum(dat$birth_time <t & dat$death_time >t)
  }),
  Death = sapply(xv, function(t){
    sum(dat$death_time <t)
  })
) |> 
  pivot_longer(Birth:Death)


ggplot(bdmf, aes(x=time, y=value, colour=name, linetype = name))+
  geom_line()+
  theme_bw()+
  scale_color_brewer(palette = "Set1")

ggsave("birth_death_marginal.png", width = 7, height = 7)

####
#MLE
ll <- function(par, ti, Tau){
  mu <- exp(par[1])
  b <- exp(par[2])
  a <- exp(par[3])
  sum(log(mu+sapply(ti, function(t)sum(a*exp(-b*(t-ti[ti<t])))))) - 
    (mu*Tau+sum(-a*expm1(-b*(Tau-ti))/b))
}

opt <- optim(rep(0,3), ll,
             ti = dat$birth_time, Tau = Tau,
             control = list(fnscale=-1), method = "BFGS")
print(opt)

muhat <- exp(opt$par[1])
bhat <- exp(opt$par[2])
ahat <- exp(opt$par[3])
print(c(muhat, ahat, bhat))

ti <- dat$birth_time

xv <- seq(0,Tau,by=0.1)
cumint_hawkes <- function(xv, mu, a, b, ti){
  sapply(xv, function(t){mu*t+sum(-a*expm1(-b*(t-ti[ti<t]))/b)})  
}

Lv <- cumint_hawkes(xv,mu,a,b,ti)
Lvhat <- cumint_hawkes(xv,muhat,ahat,bhat,ti)

dat <- mutate(dat, count = row_number())
dat_step <-  rbind(
  data.frame(id =NA, parent=NA, birth_time=0,death_time=NA, generation=NA, count=0),
  dat,
  data.frame(id =NA, parent=NA, birth_time=Tau,death_time=NA, generation=NA, count=nrow(dat)))

p1 <- ggplot() +
  geom_step(data = dat_step, aes(x=birth_time, y = count))+
  geom_line(data = NULL, aes(x=xv, y=Lv, colour="true", linetype = "true"))+
  geom_line(data = NULL, aes(x=xv, y=Lvhat, colour="MLE", linetype = "MLE"))+
  theme_classic()+
  labs(x = "time", y="cumulative couunt", linetype="param.", colour="param.")


int_hawkes <- function(xv,mu,a,b,ti){
  mu+sapply(xv, function(t)sum(a*exp(-b*(t-ti[ti<t]))))
}
lv <- int_hawkes(xv,mu,a,b,ti)
lvhat <- int_hawkes(xv,muhat,ahat,bhat,ti)

p2 <- ggplot() +
  geom_rug(data = dat, aes(x=birth_time))+
  geom_line(data = NULL, aes(x=xv, y=lv, colour="true", linetype = "true"))+
  geom_line(data = NULL, aes(x=xv, y=lvhat, colour="MLE", linetype = "MLE"))+
  theme_classic()+
  labs(x = "time", y="couunt", linetype="param.", colour="param.")

p1/p2
ggsave("intensity.png", width = 7, height = 7)

## dist. of MLE
muhat <- numeric(1000)
ahat <- numeric(1000)
bhat <- numeric(1000)
set.seed(20260820)
system.time({
  for(i in 1:1000){
    dat <- simulate_branching_death(mu, a, b, Tau)
    if(nrow(dat)>1){
      opt <- optim(rep(0,3), ll,
                   ti = dat$birth_time, Tau = Tau,
                   control = list(fnscale=-1), method = "BFGS")
      muhat[i] <- exp(opt$par[1])
      bhat[i] <- exp(opt$par[2])
      ahat[i] <- exp(opt$par[3])
    }else{
      muhat[i] <- NA_real_
      bhat[i] <- NA_real_
      ahat[i] <- NA_real_
    }
  }
})

h1 <- ggplot(data = NULL)+
  geom_histogram(aes(x=muhat - mu), fill = "lightgrey", colour="black", bins = 30)+
  theme_classic()

h2 <- ggplot(data = NULL)+
  geom_histogram(aes(x=ahat - a), fill = "lightgrey", colour="black", bins = 30)+
  theme_classic()

h3 <- ggplot(data = NULL)+
  geom_histogram(aes(x=bhat - b), fill = "lightgrey", colour="black", bins = 30)+
  theme_classic()


h1/h2/h3
ggsave("hist_mle.png", width = 7, height = 7)

####
#branching ratio

set.seed(20260820)
sim2 <- simulate_branching_death(1, 1.2, 1, 10)

segments <- dplyr::filter(sim2, parent > 0) |> 
  rowwise() |> 
  mutate(
    x = sim2$birth_time[sim2$id == parent],
    y = sim2$generation[sim2$id == parent],
    xend = birth_time,
    yend = generation-0.05
  )

ggplot(sim2) +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.02, "npc")),
    colour="grey") +
  geom_point(aes(birth_time, generation), size=3, shape=1) +
  geom_rug(aes(x=birth_time))+
  scale_y_reverse() +
  labs(x = "time", y="generation") +
  theme_bw()

ggsave("birth_death2.png", width = 7, height = 7)

bdmf <- data.frame(
  time = xv,
  Birth = sapply(xv, function(t){
    sum(sim2$birth_time <t)
  }),
  Current = sapply(xv, function(t){
    sum(sim2$birth_time <t & sim2$death_time >t)
  }),
  Death = sapply(xv, function(t){
    sum(sim2$death_time <t)
  })
) |> 
  pivot_longer(Birth:Death)


ggplot(bdmf, aes(x=time, y=value, colour=name, linetype = name))+
  geom_line()+
  theme_bw()+
  scale_color_brewer(palette = "Set1")

ggsave("birth_death_marginal2.png", width = 7, height = 7)
