library(ggplot2)

#中央値
mu <- qexp(0.5)
print(qnorm(0.5,mu,1))

#分布の形を図で確認
curve(dnorm(x,mu,1),  xlim=c(-2,6), lwd=1.5, ylab="density")
curve(dexp(x), xlim=c(0,6), add=TRUE, lty=2, col="royalblue", lwd=1.5)
legend("topright", c("normal","exp"), lty=1:2, col=c("black", "royalblue"), lwd=1.5)

#差の中央値（数値積分による力技で確認）
conv <- function(x){
  #たたみ込み
  integrate(function(y){dnorm(x+y,mu)*dexp(y)},-Inf,0)$value + 
    integrate(function(y){dnorm(x+y,mu)*dexp(y)},0,Inf)$value
}
conv <- Vectorize(conv)
m_root <- uniroot(function(x)integrate(conv,-Inf,x)$value-0.5, c(-1,0))

print(m_root$root)
#-0.1825963

#分散（数値積分による力技で確認）
v_e <- integrate(function(x)x^2*dexp(x), 0,Inf)$value - 
  integrate(function(x)x*dexp(x), 0,Inf)$value^2

v_n <- integrate(function(x)x^2*dnorm(x,mu,1), -Inf,Inf)$value -
  integrate(function(x)x*dnorm(x,mu,1), -Inf,Inf)$value^2

print(v_e)
#1
print(v_n)
#1

?wilcox.test
#繰り返し計算
n1 <-100
n2 <-100
pv_simfun <- function(n1,n2, mu, m){
  x1 <- rnorm(n1,mu,1)
  x2 <- rexp(n2)
  res_w <- wilcox.test(x1, x2, conf.int=TRUE)
  cp <- res_w$conf.int[1] < m & m < res_w$conf.int[2]
  res_w_shift <- wilcox.test(x1, x2, mu=m)
  data.frame(estimate=res_w$estimate,
             pv = res_w$p.value,
             pv2 = res_w_shift$p.value,
             cover = cp)
}

#手元のパソコンだと10秒くらいかかる
system.time({
  res <- lapply(1:10000, 
                function(i){set.seed(i);pv_simfun(25, 25, mu, m_root$root)})
  res <- do.call("rbind",res)
})
#  user  system elapsed 
#10.068   0.051  10.212 

#p値の経験分布
ggplot(res,aes(x=pv))+
  stat_ecdf()+
  geom_abline(slope = 1,intercept = 0, lty=2)+
  theme_bw()+labs(x="nominal", y="actual")+
  theme_bw()+labs(x="nominal", y="actual", title="p-value (scenario 1)")
ggsave("wilcox_p.png", width = 8, height=8)

print(round(mean(res$cover),2))
#0.95

ggplot(res,aes(x=pv2))+
  stat_ecdf()+
  geom_abline(slope = 1,intercept = 0, lty=2)+
  theme_bw()+labs(x="nominal", y="actual", title="p-value (scenario 2)")
ggsave("wilcox_p2.png", width = 8, height=8)
