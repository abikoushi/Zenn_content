library(ggplot2)
library(dplyr)

sample2 <- function(x, size) {
  if(length(x)==1) {
    rep(x, size)
  } else{
    sample(x, size, replace=TRUE)
  }
}


sir.full <- function(size,
                     R0=2.5,
                     meanlog=1.62,
                     sdlog=0.42,
                     meanlog2=1.54,
                     sdlog2=0.37,
                     I0,
                     seed = NULL,
                     imax,
                     keep.intrinsic=FALSE){
  
  if (!is.null(seed)) set.seed(seed)
  
  V <- 1:size
  
  if (missing(I0)) {
    stop("specify the initial conditions")
  }
  
  initial_infected <- 1:I0
  
  if (missing(imax)) imax <- size
  
  # queue_v : 次の感染者のID
  # queue_t : 次の感染者の感染時刻
  # queue_infector : 次の感染の感染源のID
  queue_v <- queue_t <- queue_infector <- rep(NA, I0)
  
  queue_v[1:I0] <- initial_infected 
  queue_t[1:I0] <- 0
  
  t_infected <- t_symptomatic <- rep(NA, size)
  t_infected[initial_infected] <- 0
  
  t_gillespie <- NULL
  c_infected <- 0
  
  if (keep.intrinsic) {
    intrinsic_generation <- vector('list', length(V))
  } else {
    intrinsic_generation <- NULL
  }
  
  done <- rep(FALSE, size)
  infected_by <- rep(NA, size)
  
  stop <- FALSE
  
  while (!stop) {
    j.index <- which.min(queue_t) #キューの中から感染時刻が最も早い個体を選択
    j <- queue_v[j.index]
    
    infected_by[j] <- queue_infector[j.index]
    t_infected[j] <- queue_t[j.index]
    
    t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
    
    incubation <- rlnorm(1, meanlog, sdlog)
    t_symptomatic[j] <- t+incubation #発症時刻
    
    c_infected <- c_infected + 1 #累計感染者数
    
    ncontact <- rpois(1, R0) #平均R0人の感染者が発生
    
    n <- V[V != j]
    
    if (ncontact > 0) {
      queue_v <- c(queue_v, sample2(n, ncontact)) #感染者j以外の個体から、二次感染者を抽出
      queue_infector <- c(queue_infector, rep(j, ncontact))
    }
    
    generation <- rlnorm(ncontact, meanlog2, sdlog2) #世代時間
    
    if (keep.intrinsic) intrinsic_generation[[j]] <- generation
    
    if (ncontact > 0) {
      queue_t <- c(queue_t, t + generation)
    }
    
    done[j] <- TRUE #感染済みのフラグ
    
    filter2 <- !done[queue_v] #未感染の個体をキューに残す
    queue_v <- queue_v[filter2]
    queue_infector <- queue_infector[filter2]
    queue_t <- queue_t[filter2]
    
    #終了条件
    stop <- (c_infected == length(V) || #集団全体が感染
               length(queue_v) == 0L || #キュー内に未感染の個体が存在しない
               c_infected == imax) #感染者数が imax に到達
  }
  
  return(
    list(
      data=data.frame(
        time=t_gillespie[(I0):c_infected], #感染イベントの発生時刻
        infected=(I0):c_infected #累計感染者数
      ),
      intrinsic_generation=intrinsic_generation, #世代時間
      t_infected=t_infected, #感染時刻
      t_symptomatic=t_symptomatic, #発症時刻
      infected_by=infected_by #感染源のID
    )
  )
}


## オイラー・ロトカ方程式を指数関数的増加率 r について解く
get_growth_rate <- function(R0, meanlog2, sdlog2) {
  
  f <- function(r) {
    integrate(
      function(g){exp(-r * g) * dlnorm(g, meanlog = meanlog2, sdlog = sdlog2)},
      lower = 0,
      upper = Inf
    )$value - 1 / R0
  }
  
  upper <- 1
  while (f(upper) > 0) {
    upper <- upper * 2
  }
  
  uniroot(f, interval = c(0, upper))$root
}



size <- 1000000
R0 <- 2.5

meanlog <- 1.62
sdlog <- 0.42

meanlog2 <- 1.54
sdlog2 <- 0.37

system.time({
  sim <- sir.full(
    size = size,
    R0 = R0,
    meanlog = meanlog,
    sdlog = sdlog,
    meanlog2 = meanlog2,
    sdlog2 = sdlog2,
    I0 = 1,
    imax = 10000,
    seed = 101,
    keep.intrinsic = FALSE
  )
})
#   user  system elapsed 
# 37.706   0.176  37.986 

r <- get_growth_rate(
  R0 = R0,
  meanlog2 = meanlog2,
  sdlog2 = sdlog2
)

ggplot(sim$data, aes(time,infected))+
  geom_line()+
  stat_function(fun=function(x){exp(r*x)}, colour="orangered", linetype=2)+
  ylim(c(0,11000))+
  theme_classic(12)

#ggsave(filename = "infected.png", width = 5, height = 5)

dat_pairs <- data.frame(t_symptomatic = sim$t_symptomatic,
                        t_infected = sim$t_infected,
                        infected_by = sim$infected_by) |> 
  dplyr::filter(!is.na(infected_by) & !is.na(t_infected)) |>  #感染した個体だけをフィルタ
  mutate(parent_symptomatic = sim$t_symptomatic[infected_by],
         parent_infected = sim$t_infected[infected_by]) |> 
  mutate(serial = t_symptomatic - parent_symptomatic, 
         incubation = t_symptomatic - t_infected,
         parent_incubation = parent_symptomatic - parent_infected,
         generation = t_infected - parent_infected)

print(head(dat_pairs))

dat_pairs_f <- dplyr::filter(dat_pairs, 39 <= parent_symptomatic & parent_symptomatic <= 40)
dat_pairs_b <- dplyr::filter(dat_pairs, 39 <= t_symptomatic & t_symptomatic <= 40)

####
int <- integrate(function(x){dlnorm(x, meanlog, sdlog)*exp(-r*x)},0,Inf)
p_biased <- function(x){
  sapply(x, function(xv)integrate(function(t){dlnorm(t, meanlog, sdlog)*exp(-r*t)/int$value},0,xv)$value)
}

ggplot(dat_pairs_f, aes(x=parent_incubation))+
  stat_function(fun = function(x)plnorm(x, meanlog, sdlog))+
  stat_function(fun = p_biased, colour="darkgrey")+
  stat_ecdf(aes(colour="forward", linetype = "forward")) + 
  stat_ecdf(data = dat_pairs_b, aes(colour="backward", linetype = "backward"))+
  labs(colour="method", linetype="method")+
  theme_classic(12)

ggsave(filename = "inc_cdf.png", width = 5, height = 5)

####
cat("forward:", 1/mean(exp(-r*dat_pairs_f$serial)) , "\n")
cat("backward:", 1/mean(exp(-r*dat_pairs_b$serial)) , "\n")
# > cat("forward:", 1/mean(exp(-r*dat_pairs_f$serial)) , "\n")
# forward: 2.434675 
# > cat("backward:", 1/mean(exp(-r*dat_pairs_b$serial)) , "\n")
# backward: 1.580067 

#####
#以下 Park et al (2021) Fig.2 の図（記事本文では不使用） 
#####
library(patchwork)
theme_set(theme_classic(12))

## incubation
meanlog1 <- 1.62
sdlog1 <- 0.42
# mean(rlnorm(50000, meanlog = 1.62, sdlog=0.42)) ## 5.5 days
# sd(rlnorm(50000, meanlog = 1.62, sdlog=0.42)) ## 2.4
#exp(meanlog1+sdlog1^2/2)
#sqrt( exp(2*meanlog1+sdlog1^2)*expm1(sdlog1^2) )

serialR <- function(r,
                    rho,
                    meanlog1=1.62,
                    sdlog1=0.42,
                    meanlog2=1.54,
                    sdlog2=0.37,
                    nsim=500000,
                    seed=101,
                    resampling=TRUE) {
  set.seed(seed)
  backinc <- rlnorm(nsim, meanlog=meanlog1, sdlog=sdlog1)
  if(resampling){
    backinc <- sample(backinc, replace=TRUE, prob=exp(-r*backinc))    
  }
  
  muloggen <- meanlog2 + sdlog2 * rho * (log(backinc) - meanlog1)/sdlog1
  sigmaloggen <- sqrt(sdlog2^2 * (1 - rho^2))
  
  forwardgen <- rlnorm(nsim, muloggen, sigmaloggen)
  forwardinc <- rlnorm(nsim, meanlog=meanlog1, sdlog=sdlog1)
  
  ser <- -backinc + forwardgen + forwardinc
  
  
  data.frame(
    r = r,
    rho = rho,
    R =1/mean(exp(-r*ser))
  )
}


generationR_int <- function(r,
                            meanlog2=1.54,
                            sdlog2=0.37) {
  int <- integrate(function(x){exp(-r*x)*dlnorm(x, meanlog2, sdlog2)},
                   0, Inf)
  
  return(1/int$value)
}


condf <- expand.grid(r=seq(0, 0.3, by=0.01),
                     rho=c(0,0.25,0.5,0.75))
res1 <- bind_rows(mapply(serialR, condf$r, condf$rho, SIMPLIFY = FALSE))

res_true <- distinct(condf, r) |> 
  rowwise() |> 
  mutate(R = generationR_int(r),
         rho="g")

ggplot(res1,aes(x=r,y=R, group=rho,
                colour=factor(rho),
                shape = factor(rho), 
                linetype = factor(rho)))+
  geom_point()+
  geom_line()+
  geom_line(data=res_true, aes(y=R))+
  geom_point(data=res_true, aes(y=R))+
  scale_color_brewer(palette =  "Set2")+
  labs(colour="correlation", shape="correlation", linetype="correlation")

res2 <- bind_rows(mapply(serialR, condf$r, condf$rho, resampling=FALSE,
                         SIMPLIFY = FALSE))

ggplot(res2,aes(x=r,y=R, group=rho,
                colour=factor(rho),
                shape = factor(rho), 
                linetype = factor(rho)))+
  geom_point()+
  geom_line()+
  geom_line(data=res_true, aes(y=R))+
  geom_point(data=res_true, aes(y=R))+
  scale_color_brewer(palette =  "Set2")+
  labs(colour="correlation", shape="correlation", linetype="correlation")

