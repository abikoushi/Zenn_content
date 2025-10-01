library(dqrng) #デフォルトの乱数生成より速いらしい
Ti = 10
r = 0.9
k = 11
#上側確率
res = poisson.test(x = k, T = Ti, r = r, alternative = "greater")
res$p.value
ppois(k-1, lambda = Ti*r, lower.tail = FALSE)
pgamma(r, shape = k, rate = Ti, lower.tail = TRUE)

#下側確率
res = poisson.test(x = k, T = Ti, r = r, alternative = "less")
res$p.value
ppois(k, lambda = Ti*r, lower.tail = TRUE)
pgamma(r, shape = k+1, rate = Ti, lower.tail = FALSE)

#両側
res = poisson.test(x = k, T = Ti, r = r, alternative = "two.sided")
res$conf.int
c(qgamma(0.025, shape = k, rate = Ti),
  qgamma(0.975, shape = k+1, rate = Ti))

res = poisson.test(x = k, T = Ti, r = r, alternative = "two.sided")


pvfun = function(r0,k,Ti) {
  lower = pgamma(r0, shape = k, rate = Ti, lower.tail = TRUE)
  upper = pgamma(r0, shape = k+1, rate = Ti, lower.tail = FALSE)
  2*pmin(lower, upper,0.5)
}

pvfun_b = function(r0,k,Ti, a, b) {
  lower = pgamma(r0, shape = k+a, rate = Ti+b, lower.tail = TRUE)
  upper = pgamma(r0, shape = k+a, rate = Ti+b, lower.tail = FALSE)
  2*pmin(lower, upper)
}

ci_b = c(qgamma(0.025, shape = k+0.5, rate = Ti),
    qgamma(0.975, shape = k+0.5, rate = Ti))

eroorbar = function(ci, y, length=0.1,angle=90,...){
  arrows(ci[1], y, ci[2], y, code=3, angle=angle, length = length, ...)
}

png("pvfun_pois.png", width = 500, height = 500)
curve(pvfun(x,k,Ti), 0, 3, xlab="r", ylab="alpha", n=501,lwd=1.5)
eroorbar(res$conf.int, 0.05, lty=2,lwd=1.5)
curve(pvfun_b(x,k,Ti,0.5,0),add=TRUE, n=501, col="royalblue",lwd=1.5)
eroorbar(ci_b, 0.05, lty=2,lwd=1.5, col="royalblue")
legend("topright", legend = c("poisson.test","Bayesian"), lty=1, col=c("black","royalblue"), lwd=1.5)
dev.off()

####

sim_waitingtime = function(iter, k){
  x = numeric(iter)
  for(it in seq_len(iter)){
    u = 0
    j = 1L
    while(j < k){
      ep = dqrng::dqrexp(1)
      u = u + ep
      j = j+1L
    }
    x[it] = u
  }
  return(x)
}

sim_count = function(iter, tau){
  x=integer(iter)
  for(it in seq_len(iter)){
    u=0
    j=0L
    while(u < tau){
      ep = dqrng::dqrexp(1)
      u = u + ep
      j = j+1L
    }
    x[it] = j-1L
  }
  return(x)
}

set.seed(20251001)
X = sim_count(100000, Ti)
#hist(X)
tab = table(X)
ran = range(X)
png("simcount.png",width = 500, height = 500)
plot(tab/sum(tab), ylab="Probability")
points(ran[1]:ran[2], dpois(ran[1]:ran[2], Ti), type="o", col="royalblue", lwd=2, lty=3)
dev.off()

Y = sim_waitingtime(100000, k)
png("simtime.png",width = 500, height = 500)
hist(Y, breaks = "scott", freq = FALSE, main="", xlab = "waiting time", border = "lightgrey")
curve(dgamma(x, k-1, 1), col="royalblue", add=TRUE, lwd=2, lty=2)
dev.off()
