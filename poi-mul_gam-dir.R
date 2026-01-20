library(ggplot2)

poi_to_mult <- function(N,lambda){
  X = matrix(0L, N, length(lambda))
  i = 0L
  while (i<=N) {
    x = rpois(length(lambda),lambda)
    if(sum(x) == 10L){
      X[i,] = x
      i = i+1L
    }
  }  
  return(X)
}

lambda = c(1,2,4)
system.time({
  set.seed(1234);res_mult <- poi_to_mult(10000,lambda)  
})
#  user  system elapsed 
# 0.066   0.000   0.066 

# phat = colMeans(res_mult/10)
# bp = barplot(phat, names.arg = 1:4, ylim=c(0,0.5))
# points(bp, lambda/sum(lambda), type="b")

df_mult <- reshape2::melt(res_mult)

pp = ggplot(df_mult)+
  stat_ecdf(aes(x=value, group=Var2, colour=factor(Var2)))+
  stat_function(fun = pbinom, args=list(size=10, prob=lambda[1]/sum(lambda)), linetype=2)+
  stat_function(fun = pbinom, args=list(size=10, prob=lambda[2]/sum(lambda)), linetype=2)+
  stat_function(fun = pbinom, args=list(size=10, prob=lambda[3]/sum(lambda)), linetype=2)+
  labs(colour="variable", y="ecdf")+
  theme_classic(18)

print(pp)
ggsave(filename = "mult_cdf.png", plot = pp, width=7, height = 7)
####

gamm_to_dir <- function(N, alpha){
  X = matrix(0L, N, length(alpha))
  for(i in 1:10000){
    x = rgamma(length(alpha), alpha, 1)
    X[i,] = x/sum(x)    
  }
  return(X)
}

alpha = c(1,2,4)
system.time({
  set.seed(2345); res_dir <- gamm_to_dir(10000, alpha)
})
#  user  system elapsed 
# 0.009   0.000   0.009 


df_dir <- reshape2::melt(res_dir)

pp = ggplot(df_dir)+
  stat_ecdf(aes(x=value, group=Var2, colour=factor(Var2)))+
  stat_function(fun = pbeta, args=list(shape1=alpha[1], shape2 = sum(alpha[-1])), linetype=2)+
  stat_function(fun = pbeta, args=list(shape1=alpha[2], shape2 = sum(alpha[-2])), linetype=2)+
  stat_function(fun = pbeta, args=list(shape1=alpha[3], shape2 = sum(alpha[-3])), linetype=2)+
  labs(colour="variable", y="ecdf")+
  theme_classic(18)

print(pp)
ggsave(filename = "dir_cdf.png", plot = pp, width=7, height = 7)
