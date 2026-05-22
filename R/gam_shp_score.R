library(dplyr)
library(ggplot2)

gamma_shape_score_test <- function(x, alpha0) {
  n <- length(x)
  xbar <- mean(x)
  logxbar <- mean(log(x))
  
  U <- n * (logxbar - log(xbar) + log(alpha0) - digamma(alpha0))
  Ieff <- n * (trigamma(alpha0) - 1 / alpha0)
  stat <- (U^2) / Ieff
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  
  list(
    statistic = stat,
    p.value = pval,
    null.value = alpha0
  )
}

pv_simfun <- function(n, shape0, shape, scale){
  pv <- numeric(10000)
  for(i in 1:10000){
    x <- rgamma(n, shape = shape, scale = scale)
    pv[i] <- gamma_shape_score_test(x, shape0)$p.value  
  }
  return(pv)
}

set.seed(987654321)
conddf <- expand.grid(n=c(10, 25, 50, 100), 
                      shape=c(0.5, 1, 2),
                      scale=1)
res_alpha <- rowwise(conddf) %>% 
  reframe(n=n, shape=shape, scale=scale, 
          pv=pv_simfun(n, shape, shape, scale))

conddf <- expand.grid(n=c(10, 25, 50, 100), 
                      shape=c(0.5, 0.75, 1.5, 2),
                      scale=1)

res_beta <- rowwise(conddf) %>% 
  reframe(n=n, shape=shape, scale=scale, 
          pv=pv_simfun(n, 1, shape, scale))


ggplot(res_alpha, aes(x=pv))+
  geom_abline(slope = 1, intercept=0, linetype=2, colour="royalblue")+
  stat_ecdf()+
  facet_grid(n~shape, labeller = label_both)+
  theme_bw(16)+labs(x="nominal", y="actual")
ggsave("alpha.png", width = 10, height = 7)

ggplot(res_beta, aes(x=pv))+
  geom_abline(slope = 1, intercept=0, linetype=2, colour="royalblue")+
  stat_ecdf()+
  facet_grid(n~shape, labeller = label_both)+
  theme_bw(16)+labs(x="nominal", y="actual")
ggsave("beta.png", width = 12, height = 7)
