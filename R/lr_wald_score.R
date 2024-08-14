#尤度比検定
LRtest <- function(p,x,n){
  phat <- x/n
  lp <- 2*(dbinom(x, n, phat ,log = TRUE)-
             dbinom(x, n, p ,log = TRUE))
  pchisq(lp, 1, lower.tail = FALSE)
}
set.seed(123)
n <- 20
x <- rbinom(1, n, 0.5)
print(x)
# [1] 9

z <- seq(0.01, 0.99, by=0.005)
pv_lr <- sapply(z, LRtest, x=x, n=n)
df_lr <- data.frame(p=z, pv=pv_lr, method="LR")

library(ggplot2)
ggplot(df_lr, aes(x=p, y=pv))+
  geom_line()+
  theme_bw(18)+labs(x="param.", y="p-value")

ggsave("pfun_lr0.png")

library(rootSolve)
sol <- uniroot.all(function(p)LRtest(p,x,n)-0.05, c(0.1,0.8))
print(sol)
ggplot(df_lr, aes(x=p, y=pv))+
  geom_line()+
  geom_errorbarh(data = NULL, aes(xmin=sol[1], xmax=sol[2], y=0.05), height=0.03, colour="cornflowerblue")+
  theme_bw(18)+labs(x="param.", y="p-value")

ggsave("pfun_lr.png")

##
waldtest <- function(p,x,n){
  phat <- x/n
  se <- sqrt(phat*(1-phat)/n)
  2*pnorm(abs(phat-p)/se, lower.tail = FALSE)
}

pv_w <- sapply(z, waldtest, x=x, n=n)
df_w <- data.frame(p=z, pv=pv_w, method="Wald")

####
#Score test

res_p <- prop.test(x, n, p = 0.5, correct = FALSE, conf.level = 0.95)

CI_score <- function(x, n, level){
  z <- qnorm(0.5*(1-level), lower.tail = FALSE)
  phat <- x/n
  t_1 <- phat+(z^2)/(2*n)
  t_2 <- z*sqrt(z^2/(4*n^2)+phat*(1-phat)/n)
  c((t_1 - t_2)/(1+(z^2)/n),
    (t_1 + t_2)/(1+(z^2)/n))  
}


print(CI_score(x,n,0.95))
print(res_p$conf.int)

pv_p <- sapply(z, function(p)prop.test(x, n, p = p, correct = FALSE)$p.value)
df_p <- data.frame(p=z, pv=pv_p, method="score")
df_pv <- rbind(df_lr,df_w,df_p)
ggplot(data = df_pv, aes(x=p, y=pv, colour=method, group = method, linetype=method))+
  geom_line()+
  scale_color_brewer(palette = "Set2")+
  theme_bw(16)+
  labs(x="param.", y="p-value", colour="method", linetype="method")
ggsave("proptest.png")