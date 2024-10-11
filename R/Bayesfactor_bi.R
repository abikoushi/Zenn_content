# x <- 7
# n <- 10
# #check the numerical solution
# integrate(function(p){(p^x)*((1-p)^(n-x))}, lower = 0, upper = 1)
# #check the analytic solution
# beta(x+1, n-x+1)

BF_bi <- function(x, n, p0=0.5){
   (x*log(p0)+(n-x)*log(p0))-lbeta(x+1, n-x+1)
}
n <- 10
x <- rbinom(10000, size = n, p=0.5)
bstat <- BF_bi(x,n)
range(bstat)
hist(bstat)
