library(ggplot2)
library(patchwork)
library(rootSolve)
library(MCMCpack)

pb <- dbinom(0:10, 10, 0.5)
t_p <- dbinom(8, 10, 0.5)
cols = ifelse(pb <= t_p, "orangered", "grey80")
#png("binom_barplot.png")
barplot(pb, names.arg = 0:10, col=cols)
#dev.off()
cat("p-value (by definition): ", sum(pb[pb <= t_p]), "\n")

res0_b <- binom.test(x = 8, n = 10, p = 0.5)
cat("p-value (binom.test): ", res0_b$p.value, "\n")
print(res0_b$p.value)

z <- seq(0.1, 0.99, by=0.005)
pv_b <- sapply(z, function(p)binom.test(8, 10, p = p)$p.value)

ggplot(data = NULL)+
  geom_line(aes(x=z, y=pv_b))+
  geom_errorbarh(aes(y=0.05, xmin = res0_b$conf.int[1], xmax=res0_b$conf.int[2]),
                 height=0.03, colour="cornflowerblue")+
  theme_bw(16)+
  labs(x="param.", y="p-value", colour="method", linetype="method")
#ggsave("pvfun_binom.png")

###
#binom.test
pvfun0 <- function(p, x, n){
  pu <- pbinom(x-1, n, p, lower.tail = FALSE)
  pl <- pbinom(x, n, p, lower.tail = TRUE)
  2*pmin(pl, pu)
}
pv_b0 <- pvfun0(z, 8, 10)
ggplot(data = NULL)+
  geom_line(aes(x=z, y=pv_b, colour="by definicion"))+
  geom_line(aes(x=z, y=pv_b0, colour="twice one-side"))+
  geom_errorbarh(aes(y=0.05, xmin = res0_b$conf.int[1], xmax=res0_b$conf.int[2], colour="twice one-side"),
                 height=0.03)+
  theme_bw(16)+
  labs(x="param.", y="p-value", colour="method", linetype="method")
ggsave("pvfun_binom2.png")

###
#poisson.test
pvfun0 <- function(r, x, tau){
  pu <- ppois(x-1, r*tau, lower.tail = FALSE)
  pl <- ppois(x, r*tau, lower.tail = TRUE)
  2*pmin(pl, pu)
}
rvec <- seq(0.01,3,by=0.01)
pv_p <- sapply(rvec, function(r)poisson.test(8, 10, r = r)$p.value)
pv_p0 <- pvfun0(rvec, 8, 10)
res0_p <- poisson.test(8, 10)

ggplot(data = NULL)+
  geom_line(aes(x=rvec, y=pv_p, colour="by definicion"))+
  geom_line(aes(x=rvec, y=pv_p0, colour="twice one-side"))+
  geom_errorbarh(aes(y=0.05, xmin = res0_p$conf.int[1], xmax=res0_p$conf.int[2], colour="twice one-side"),
                 height=0.03)+
  theme_bw(16)+
  labs(x="param.", y="p-value", colour="method", linetype="method")

ggsave("pvfun_pois.png")

#####
#fisher.test
confint_fisher <- function(X, level=0.95){
  rs <- rowSums(X)
  cs <- colSums(X)
  alpha <- 1-level
  testfun_up <- function(phi){
    pall <- MCMCpack::dnoncenhypergeom(x = NA, cs[1],cs[2],rs[1], exp(phi))
    i <- match(X[1,1],pall[,1])
    pv <- sum(pall[i:nrow(pall),2])
    pv-alpha/2
  }
  testfun_low <- function(phi){
    pall <- MCMCpack::dnoncenhypergeom(x = NA,
                                       cs[1],cs[2],rs[1], exp(phi))
    i <- match(X[1,1],pall[,1])
    pv <- sum(pall[1:i,2])
    pv-alpha/2
  }
  res_u <- uniroot(testfun_up, lower = -10, upper = 10)
  res_l <- uniroot(testfun_low, lower = -10, upper = 10)
  list(oddsratio=c(exp(res_u$root),exp(res_l$root)))
}
X <-matrix(c(12,5,6,12), nrow=2)

resf <- fisher.test(X, conf.level = 0.95)
print(resf$conf.int)
print(confint_fisher(X))

exact2x2::exact2x2(X)

###
#まとめ

pvfun <- Vectorize(FUN = function(mu0){
  binom.test(8, 10,  p=mu0)$p.value  
})
sol_ci <- rootSolve::uniroot.all(f = function(p){pvfun(p)-0.05}, interval = c(0.1,0.99))
print(sol_ci)

png("pfun_simple.png")
curve(pvfun(x), 0.1, 0.99, n=1001)
dev.off()

##


