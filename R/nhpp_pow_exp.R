NHPP_exp <- function(n, a, b) {
  mlogz <- -log(runif(n))
  z <- exp(b) + a * cumsum(mlogz)
  t <- rep(Inf, n)
  idx <- z > 0
  t[idx] <- (log(z[idx]) - b) / a
  t
}

Weibull_process <- function(n, alpha, beta) {
  z <- runif(n)
  alpha * cumsum(-log(z))^(1 / beta)
}
