library(rootSolve)
library(exact2x2)

pvfun0 <- function(p, x, n){
  pu <- pbinom(x-1, n, p, lower.tail = FALSE)
  pl <- pbinom(x, n, p, lower.tail = TRUE)
  2*pmin(pl, pu)
}
sol0 <-uniroot.all(f = function(p){pvfun0(p,5,10)-0.05},
            interval = c(0.1,0.9))
res_b <- binom.test(5, 10)
print(res_b$conf.int)
print(sol0)

pvfun <- function(p){
  sapply(p, function(p)binom.test(5,10,p)$p.value)
}
sol <- uniroot.all(f = function(p){pvfun(p)-0.05},
            interval = c(0.1,0.9))
print(sol)

curve(pvfun0(x, 5, 10))
curve(pvfun(x), col="red", add=TRUE)

###
pvec <- seq(0.01,0.99,by=0.01)
pval <- sapply(pvec, function(p)binom.test(5, 10, p = p)$p.value)
#png("confint1.png")
plot(pvec, pval, type="s", xlab="param", ylab="p-value")

segments(x0 = res_b$conf.int[1], x1 = res_b$conf.int[2],
         y0 = 0.05, y1 = 0.05, col="royalblue")
#dev.off()
png("confint2.png")
plot(pvec, pval, type="s", xlab="param", ylab="p-value")
segments(x0 = sol[1], x1 = sol[2],
         y0 = 0.05, y1 = 0.05, col="orangered")
dev.off()


###
library(rootSolve)
pvfun0 <- function(r, x, tau){
  pu <- ppois(x-1, r*tau, lower.tail = FALSE)
  pl <- ppois(x, r*tau, lower.tail = TRUE)
  2*pmin(pl, pu)
}
sol0 <-uniroot.all(f = function(p){pvfun0(p,5,10)-0.05},
                   interval = c(0.1,3))
sol0
res_p
res_p <- poisson.test(5, 10, r=1)
print(res_p$conf.int)
print(sol0)

#confidence interval from p-value 
pvfun <- function(r){
  sapply(r, function(r)poisson.test(5,10,r=r)$p.value)
}
sol <- uniroot.all(f = function(p){pvfun(p)-0.05},
                   interval = c(0.1,3))
print(sol)
###
#plot p-value fun
rvec <- seq(0.01,3,by=0.01)
pval <- sapply(rvec, function(r)poisson.test(5, 10, r = r)$p.value)
png("confint1.png")
plot(rvec, pval, type="s", xlab="param", ylab="p-value")
segments(x0 = res_p$conf.int[1], x1 = res_p$conf.int[2],
         y0 = 0.05, y1 = 0.05, col="royalblue")
dev.off()
png("confint2.png")
plot(rvec, pval, type="s", xlab="param", ylab="p-value")
segments(x0 = sol[1], x1 = sol[2],
         y0 = 0.05, y1 = 0.05, col="orangered")
dev.off()

###

pfun <- Vectorize(FUN = function(mu0){
  dat <- c(2,5,0,1,3)
  t.test(dat, mu=mu0)$p.value  
})
png("pfun.png")
curve(pfun(x), -2, to = 6, n=1001)
dev.off()
