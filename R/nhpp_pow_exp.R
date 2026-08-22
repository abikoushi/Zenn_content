library(dplyr)
library(ggplot2)
library(patchwork)

NHPP_power <- function(n, alpha, beta) {
  z <- runif(n)
  alpha * cumsum(-log(z))^(1 / beta)
}

Intensity_power <- function(x,alpha,beta){
  (x/alpha)^beta
}

set.seed(123456)
res1 <- lapply(1:100, function(i)data.frame(id=i, x = NHPP_power(100, 1, 2), y=1:100))
res1 <- bind_rows(res1)
res2 <- lapply(1:100, function(i)data.frame(id=i, x = NHPP_power(100, 1, 0.5), y=1:100))
res2 <- bind_rows(res2)

p1 <- ggplot(res1, aes(x=x, y=y))+
  geom_step(aes(group=id), alpha=0.1)+
  stat_function(fun = Intensity_power, args = list(alpha=1, beta=2),
                colour="royalblue", 
                linewidth=1,
                linetype=2)+
  ylim(c(0,100))+
  theme_bw()+ggtitle(label = "beta=2")

p2 <- ggplot(res2, aes(x=x, y=y))+
  geom_step(aes(group=id), alpha=0.1)+
  stat_function(fun = Intensity_power, args = list(alpha=1, beta=0.5),
                colour="royalblue", 
                linewidth=1,
                linetype=2)+
  ylim(c(0,100))+
  theme_bw()+ggtitle(label = "beta=0.5")

p <- p1+p2
print(p)
ggsave("nhpp_pow.png", width = 7, height = 7)

#####

NHPP_exp <- function(n, a, b) {
  mlogz <- -log(runif(n))
  z <- exp(b) + a * cumsum(mlogz)
  t <- rep(Inf, n)
  idx <- z > 0
  t[idx] <- (log(z[idx]) - b) / a
  t
}

Intensity_exp <- function(t,a,b){
  (exp(a*t+b)-exp(b))/a
}

set.seed(123456)
res1 <- lapply(1:100, function(i)data.frame(id=i, x = NHPP_exp(100, 1, 0.1), y=1:100))
res1 <- bind_rows(res1)
res2 <- lapply(1:100, function(i)data.frame(id=i, x = NHPP_exp(100, -0.01, 0.1), y=1:100))
res2 <- bind_rows(res2)

p1 <- ggplot(res1, aes(x=x, y=y))+
  geom_step(aes(group=id), alpha=0.1)+
  stat_function(fun = Intensity_exp, args = list(a=1, b=0.1), 
                colour="royalblue",
                linewidth=1,
                linetype=2)+
  ylim(c(0,100))+
  theme_bw()+ggtitle(label = "a=1")

p2 <- ggplot(res2, aes(x=x, y=y))+
  geom_step(aes(group=id), alpha=0.1)+
  stat_function(fun = Intensity_exp, args = list(a=-0.01, b=0.1),
                colour="royalblue", 
                linewidth=1,
                linetype=2)+
  ylim(c(0,100))+
  theme_bw()+ggtitle(label = "a=-0.01")

p <- p1+p2
print(p)

ggsave("nhpp_exp.png", width = 7, height = 7)

