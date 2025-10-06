library(dqrng)
n = 10L
k = 4L
p = 0.6

# "greater"
pbinom(k-1, n, p, lower.tail = FALSE)
pbeta(p, shape1 = k, shape2 = n-k+1)
binom.test(k, n, p=p, alternative = "greater")$p.value

# "less"
pbinom(k, n, p)
pbeta(p, shape1 = k+1, shape2 = n-k, lower.tail = FALSE)
binom.test(k, n, p=p, alternative = "less")$p.value

res = binom.test(k, n, p=p)
res$conf.int
c(qbeta(0.025, shape1 = k, shape2 = n-k+1),
  qbeta(0.975, shape1 = k+1, shape2 = n-k))

ci_b = c(qbeta(0.025, shape1 = k+1, shape2 = n-k+1),
  qbeta(0.975, shape1 = k+1, shape2 = n-k+1))

eroorbar = function(ci, y, length=0.1,angle=90,...){
  arrows(ci[1], y, ci[2], y, code=3, angle=angle, length = length, ...)
}

pvfun = function(p0,k,N) {
  lower = pbeta(p0, shape1 = k, shape2 = N-k+1, lower.tail = TRUE)
  upper = pbeta(p0, shape1 = k+1, shape2 = N-k, lower.tail = FALSE)
  2*pmin(lower, upper,0.5)
}

pvfun_b = function(p0, k, N, a, b) {
  lower = pbeta(p0, shape1 = k+a, shape2 = N-k+b, lower.tail = TRUE)
  upper = pbeta(p0, shape1 = k+a, shape2 = N-k+b, lower.tail = FALSE)
  2*pmin(lower, upper)
}
n = 10L
k = 4L
p = 0.6
res = binom.test(k, n, p=p)

ci_b = c(qbeta(0.025, shape1 = k+1, shape2 = n-k+1),
         qbeta(0.975, shape1 = k+1, shape2 = n-k+1))

png("pvfun_binom1.png", width = 500, height = 500)
curve(pvfun(x,k,n), 0, 1, xlab="r", ylab="alpha", n=501,lwd=1.5)
eroorbar(res$conf.int, 0.05, lty=2,lwd=1.5)
curve(pvfun_b(x,k,n,1,1),add=TRUE, n=501, col="royalblue",lwd=1.5)
eroorbar(ci_b, 0.05, lty=2,lwd=1.5, col="royalblue")
legend("topright", legend = c("poisson.test","Bayesian"), lty=1, col=c("black","royalblue"), lwd=1.5)
dev.off()

n = 100L
k = 40L
res = binom.test(k, n, p = p)
ci_b = c(qbeta(0.025, shape1 = k+1, shape2 = n-k+1),
         qbeta(0.975, shape1 = k+1, shape2 = n-k+1))

png("pvfun_binom2.png", width = 500, height = 500)
curve(pvfun(x,k,n), 0, 1, xlab="r", ylab="alpha", n=1001,lwd=1.5)
eroorbar(res$conf.int, 0.05, lty=2,lwd=1.5)
curve(pvfun_b(x,k,n,1,1),add=TRUE, n=1001, col="royalblue",lwd=1.5)
eroorbar(ci_b, 0.05, lty=2,lwd=1.5, col="royalblue")
legend("topright", legend = c("poisson.test","Bayesian"), lty=1, col=c("black","royalblue"), lwd=1.5)
dev.off()
####

####
set.seed(1234)
x = runif(10)
stripchart(x)

simorderstat_unif <- function(iter, k, p, n){
  res_p = numeric(iter)
  res_k = integer(iter)
  for(it in seq_len(iter)){
    u = dqrng::dqrunif(n) #faster than stats::runif
    res_p[it] = sort(u)[k]
    res_k[it] = sum(u<p)
  }
  return(data.frame(waitingtime=res_p, counts=res_k))
}


n = 10L
k = 4L
p = 0.6
set.seed(1016)
res_sim = simorderstat_unif(100000, k, p, n)
tab = table(res_sim$counts)
ran = range(res_sim$counts)
png("simcount_binom.png", width = 500, height = 500)
plot(tab/sum(tab), ylab="Probability", xlab = "counts")
points(ran[1]:ran[2], dbinom(ran[1]:ran[2], n, p), type="o", col="royalblue", lwd=2, lty=3)
dev.off()

png("simtime_beta.png",width = 500, height = 500)
hist(res_sim$waitingtime, breaks = "scott", freq = FALSE, main="", xlab = "waiting time", border = "lightgrey")
curve(dbeta(x, k, 10-k+1), col="royalblue", add=TRUE, lwd=2, lty=2)
dev.off()
