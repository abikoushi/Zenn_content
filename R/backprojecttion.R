
lambdaest <- function(count_obs, ft, maxit=1000,tol=0.1){
  len <- length(count_obs)
  lamhat <- rep(mean(count_obs),len)
  den <-rev(cumsum(ft))
  mu <- sapply(1:len, function(t)sum(lamhat[1:t]*f[t:1]))
  ll <- numeric(maxit)
  ll[1] <- sum(dpois(count_obs,mu,log = TRUE))
  for (i in 2:maxit) {
    #E stpp
    yp <- lapply(1:len, function(t)count_obs[t]*lamhat[1:t]*ft[t:1]/sum(lamhat[1:t]*ft[t:1]))
    #M step
    num <- sapply(1:len,function(t)sum(sapply(yp, function(x)x[t]),na.rm = TRUE))
    lamhat <- (num)/(den)
    mu <- sapply(1:len, function(t)sum(lamhat[1:t]*ft[t:1]))
    ll[i] <- sum(dpois(count_obs,mu,log = TRUE))
    if(ll[i]-ll[i-1]<tol){
      break
    }
  }
  return(list(lambda=lamhat,mu=mu,it=i,ll=ll[1:i]))
}


genefun_c <- function(maxrd, maxed, mined, lambda, m, eta){
  np <- rpois(1,(maxed-mined)*lambda)
  ti <- sort(runif(np,mined,maxed))
  si <- rweibull(np,m,eta)
  rd <- ceiling(ti+si)
  id <- ceiling(ti)
  infection <- sapply(1:maxrd,function(t)sum(id==t))
  obs <- sapply(1:maxrd,function(t)sum(rd==t))
  data.frame(obs,infection)
}

set.seed(123);dat <- genefun_c(40, 20, 10, 100, 2,4)

f <- diff(pweibull(0:40,2,4))

lambda_out1 <- lambdaest(dat$obs, f, maxit  = 500, tol=0.1)
lambda_out2 <- lambdaest(dat$obs, f, maxit  = 500, tol=1e-8)
plot(lambda_out1$ll,type="l", ylab="log-likelihood")
plot(lambda_out2$ll,type="l", ylab="log-likelihood")

png("fit1.png")
plot(dat$obs,ylim=c(0,max(lambda_out$mu,dat$obs)),ylab = "sympton on set", type="h")
points(lambda_out1$mu,ylim=c(0,max(lambda_out$mu,dat$obs)),ylab = "sympton on set")
dev.off()

png("est1.png")
plot(lambda_out1$lambda,ylim=c(0,max(lambda_out$lambda,dat$infection)),ylab = "infection")
points(dat$infection,type = "s",lty=2)
dev.off()

png("est2.png")
plot(lambda_out2$lambda,ylim=c(0,max(lambda_out$lambda,dat$infection)),ylab = "infection")
points(dat$infection,type = "s",lty=2)
dev.off()



png("cumulative_intensity.png", width = 900)
par(mfrow=c(1,2))
plot(cumsum(lambda_out1$lambda),ylab = "infection", type="s", main="tol=0.1")
lines(cumsum(dat$infection),type = "s",lty=2)
plot(cumsum(lambda_out2$lambda),ylab = "infection", type="s", main="tol=0.00000001")
lines(cumsum(dat$infection),type = "s",lty=2)
dev.off()
