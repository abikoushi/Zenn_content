library(survival)
df1 <- data.frame(
  t = c(0,1,2,3,4),
  D = c(0,20,50,10,10),
  R = c(100,80,30,20,10)  
)

haz <- df1$D[-1]/df1$R[-length(df1$R)]
print(cumprod(1-haz))

time = c(rep(df1$t, df1$D), rep(df1$t[5], df1$R[1] - sum(df1$D)))
d = c(rep(1, sum(df1$D)), rep(0, df1$R[1] - sum(df1$D)))

fit1 <- survfit(Surv(time,d)~1)

#plot(fit1, conf.int = FALSE)
print(summary(fit1))

####
df2 <- data.frame(
  t = c(0,1,2,3,4),
  D = c(0,20,50,10,1),
  C = c(0,5,10,0,4),
  R = c(100,75,15,5,0)  
)

time = c(rep(df2$t, df2$D), rep(df2$t, df2$C))
d = c(rep(1, sum(df2$D)), rep(0, sum(df2$C)))

fit2 <- survfit(Surv(time,d)~1)
print(summary(fit2))

haz <- df2$D[-1]/df2$R[-length(df2$R)]
print(cumprod(1-haz))
