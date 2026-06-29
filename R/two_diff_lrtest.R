library(ggplot2)
library(dplyr)
library(tidyr)

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
  data.frame(method1=pv1, method2=pv2)
}


con_alpha <- expand.grid(n=c(10,20,100), mu=0)
res_alpha <- vector("list",3)
set.seed(1111)
for(i in 1:3){
  res <- pvsimfun(con_alpha$n[i], con_alpha$mu[i], con_alpha$mu[i])
  res_alpha[[i]] <- mutate(res, n= df_alpha$n[i])
}
res_df <- bind_rows(res_alpha)
res_df <-pivot_longer(res_df, method1:method2, names_to = "method")
ggplot(res_df, aes(x=value))+
  stat_ecdf()+
  geom_abline(slope = 1, intercept = 0, linetype=2, colour="royalblue")+
  facet_grid(n~method, labeller = label_both)+
  theme_classic(14)+labs(x="P-value")
ggsave("pv_alpha.png", width = 7, height = 7)

con_beta <- expand.grid(n=c(10,20,100), mu=c(0.1,0.2,0.5), mu0=0)
res_beta <- vector("list",9)
set.seed(2222)
for(i in 1:9){
  res <- pvsimfun(con_beta$n[i], con_beta$mu[i], con_beta$mu0[i])
  res_beta[[i]] <- mutate(res, n= con_beta$n[i], mu=con_beta$mu[i])
}
res_df <- bind_rows(res_beta)
res_df <-pivot_longer(res_df, method1:method2, names_to = "method")
ggplot(res_df, aes(x=value, colour=method, linetype = method))+
  stat_ecdf()+
  facet_grid(n~mu, labeller = label_both)+
  theme_classic(14)+labs(x="P-value")
ggsave("pv_beta.png", width = 10, height = 7)
