library(ggplot2)
library(dplyr)
library(animation)

cibeta <- function(p, x, n, a=0.5, b=0.5, level=0.95){
  alpha = 0.5*(1-level)
  ahat <- x + a
  bhat <- n - x + b
  ql = qbeta(alpha, ahat, bhat, lower.tail=TRUE)
  qu = qbeta(alpha, ahat, bhat, lower.tail=FALSE)
  c(ql, qu)
}

probbeta <- function(p, x, n, a=0.5, b=0.5){
  ahat <- x + a
  bhat <- n - x + b
  pl <- pbeta(p, ahat, bhat, lower.tail=TRUE)
  pu <- pbeta(p, ahat, bhat, lower.tail=FALSE)
  2*pmin(pl,pu)
}

z = seq(0.001, 0.999, length.out=301)

pv_beta <- probbeta(z,6,10)
CI <- cibeta(z,6,10)
df_beta <- data.frame(p=z, pv=pv_beta)

p1 = ggplot(data = df_beta, 
       aes(x=p, y=pv))+
  geom_line()+
  geom_errorbarh(data = NULL, aes(xmin = CI[1], xmax = CI[2], y=0.05), height=0.03, colour="cornflowerblue")+
  labs(y = "p-value")+
  theme_bw(16)

print(p1)
ggsave("ciplot.png", plot = p1, width = 9, height = 9)
#####

LRtest <- function(p,x,n){
  phat <- x/n
  lp <- 2*(dbinom(x, n, phat ,log = TRUE)-
             dbinom(x, n, p ,log = TRUE))
  pchisq(lp, 1, lower.tail = FALSE)
}

waldtest <- function(p,x,n){
  phat <- x/n
  se <- sqrt(phat*(1-phat)/n)
  2*pnorm(abs(phat-p)/se, lower.tail = FALSE)
}


pltfun_binom = function(n){
  plist <- vector("list",n-1L)
  z = seq(0.001, 0.999, length.out=301)
  for(x in 1L:(n-1L)){
    pv_p <- sapply(z, function(p)prop.test(x, n, p = p, correct = FALSE)$p.value)
    df_p <- data.frame(p=z, pv=pv_p, method="score")
    
    pv_w <- sapply(z, waldtest, x=x, n=n)
    df_w <- data.frame(p=z, pv=pv_w, method="Wald")
    
    pv_p <- sapply(z, function(p)prop.test(x, n, p = p, correct = FALSE)$p.value)
    df_p <- data.frame(p=z, pv=pv_p, method="score")
    
    pv_b <- sapply(z, function(p)binom.test(x, n, p = p)$p.value)
    df_b <- data.frame(p=z, pv=pv_b, method="exact")
    
    pv_lr <- sapply(z, LRtest, x=x, n=n)
    df_lr <- data.frame(p=z, pv=pv_lr, method="LR")
    
    pv_beta1 <- probbeta(z,x,n)
    df_beta1 <- data.frame(p=z, pv=pv_beta1, method="Bayesian (Jeffreys prior)")
    
    pv_beta2 <- probbeta(z,x,n,1,1)
    df_beta2 <- data.frame(p=z, pv=pv_beta2, method="Bayesian (flat prior)")
    df_pv <- bind_rows(df_lr, df_w, df_p, df_b, df_beta1, df_beta2)
    
    df_pv_s <- group_by(df_pv) %>% 
      slice(seq(1L, n(), by=6L))
    
    plist[[x]] = ggplot(data = df_pv, aes(x=p, y=pv, colour=method, shape  = method, group = method))+
      geom_line()+
      geom_point(data = df_pv_s)+
      labs(y = "p-value", title = paste("n = ", n))+
      theme_bw(16)
  }
  plist <- c(plist, rev(plist))  
  return(plist)
}

plist = pltfun_binom(10)
saveGIF({
  for (i in seq_len(length(plist))) {print(plist[[i]])}
  }, "pvalfun1.gif" , interval=0.1)

plist = pltfun_binom(20)
saveGIF({
  for (i in seq_len(length(plist))) {print(plist[[i]])}
}, "pvalfun2.gif" , interval=0.1)


plist = pltfun_binom(50)
saveGIF({
  for (i in seq_len(length(plist))) {print(plist[[i]])}
}, "pvalfun3.gif" , interval=0.1)


