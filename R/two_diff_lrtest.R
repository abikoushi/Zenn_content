library(ggplot2)
library(dplyr)

pvsimfun <- function(n, mu, mu0){
  pv1 <- numeric(10000)
  pv2 <- numeric(10000)
  for(i in 1:10000){
    x <- rnorm(n, mu)
    xbar <- mean(x)
    pv1[i] <- pnorm(sqrt(n)*(xbar-mu0), lower.tail=FALSE)
    lr2 <- 2*(sum(dnorm(x, xbar, 1, log=TRUE)) - sum(dnorm(x, mu0, 1, log=TRUE)))
    pv2[i] <- pchisq(lr2, df=1, lower.tail = FALSE)
  }
  data.frame(pv1, pv2)
}


res <- pvsimfun(10, 0, 0)
res
