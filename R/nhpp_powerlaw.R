library(animation)

intensity_powerlaw <- function(t,alpha, beta){
  (beta/alpha)*(t/alpha)^(beta-1)
}

Intensity_powerlaw <- function(x,alpha,beta){
  (x/alpha)^beta
}

NHPP_powerlaw = function(Tmax, alpha, beta, maxit){
  mlogz = rexp(1)
  t <- rep(Inf, maxit)
  t[1] = alpha*mlogz^(1/beta[1])
  for(i in 2:maxit){
    mlogz = rexp(1)
    ti = alpha*( mlogz + (t[i-1]/alpha)^beta )^( 1/beta )
    if(ti > Tmax){
      break
    }
    t[i] <- ti
  }
  return(t[1:(i-1)])
}


estimator <- function(dat, t_max){
  n <-length(dat)
  beta <- n/sum(log(t_max)-log(dat))
  alpha <- t_max/(n^(1/beta))
  c("alpha"=alpha, "beta"=beta)
}

est_intensity <- function(ti, w, Tmax){
  s <- seq(0, Tmax, by=w)
  lamhat <- sapply(seq_along(s[-1]), function(i)sum(s[i] < ti & ti < s[i+1]))/w
  data.frame(t=s[-1],intensity=lamhat)
}


###
png("intensity.png", width = 700, height=500)
par(mfrow=c(1,2))
curve(intensity_powerlaw(x,1,1), 0, 5, xlab="time", ylab = "intensity", lwd=1.5, ylim=c(0,2))
curve(intensity_powerlaw(x,1,1.5),add=TRUE, col="steelblue",lty=2, lwd=1.5, n=1001)
curve(intensity_powerlaw(x,1,0.8),add=TRUE, col="firebrick",lty=3, lwd=1.5, n=1001)
legend("topright",c("beta=1","beta=1.5","beta=0.8"), 
       lty=1:3, col=c("black","steelblue","firebrick"), lwd=1.5)

curve(Intensity_powerlaw(x,1,1), 0, 5, xlab="time", ylab = "cumulative intensity", lwd=1.5)
curve(Intensity_powerlaw(x,1,1.5),add=TRUE, col="steelblue",lty=2, lwd=1.5)
curve(Intensity_powerlaw(x,1,0.8),add=TRUE, col="firebrick",lty=3, lwd=1.5)
legend("topleft",c("beta=1","beta=1.5","beta=0.8"), 
       lty=1:3, col=c("black","steelblue","firebrick"), lwd=1.5)
dev.off()

###
Tmax <- 50
a <- 1
b <- 1.5
set.seed(1234); ti <- NHPP_powerlaw(Tmax, a, b, 1000)
est_lam <- est_intensity(ti,1,Tmax)

Ts <- seq(5, Tmax, by=1)
saveGIF({
  par(mfrow=c(1,2))
  for( i in seq_along(Ts) ){
    Tmax_c <- Ts[i]
    ti_c <- ti[ti < Tmax_c]
    est <- estimator(ti_c, Tmax_c)
    plot(c(0, ti), c(0, seq_along(ti)), type="s", xlab="time", ylab = "cumulative intensity", col="grey")
    lines(c(0, ti_c), c(0, seq_along(ti_c)), type="s")
    curve(Intensity_powerlaw(x,est[1],est[2]), add=TRUE, col="orange")
    curve(Intensity_powerlaw(x,a,b), add=TRUE, col="royalblue", lty=2)
    
    est_c <- est_intensity(ti_c, 1, Tmax_c)
    plot(est_lam$t, est_lam$intensity, type="l", xlab="time", ylab = "intensity", col="grey")
    lines(est_c$t, est_c$intensity)
    curve(intensity_powerlaw(x,est[1],est[2]), add=TRUE, col="orange")
    curve(intensity_powerlaw(x,a,b), add=TRUE, col="royalblue", lty=2)
  }
}, movie.name = "animation1.gif", interval=0.1)

Tmax <- 50
a <- 1
b <- 0.8
set.seed(1234); ti <- NHPP_powerlaw(Tmax, a, b, 1000)
est_lam <- est_intensity(ti,1,Tmax)

Ts <- seq(5, Tmax, by=1)
saveGIF({
  par(mfrow=c(1,2))
  for( i in seq_along(Ts) ){
    Tmax_c <- Ts[i]
    ti_c <- ti[ti < Tmax_c]
    est <- estimator(ti_c, Tmax_c)
    plot(c(0, ti), c(0, seq_along(ti)), type="s", xlab="time", ylab = "cumulative intensity", col="grey")
    lines(c(0, ti_c), c(0, seq_along(ti_c)), type="s")
    curve(Intensity_powerlaw(x,est[1],est[2]), add=TRUE, col="orange")
    curve(Intensity_powerlaw(x,a,b), add=TRUE, col="royalblue", lty=2)
    
    est_c <- est_intensity(ti_c, 1, Tmax_c)
    plot(est_lam$t, est_lam$intensity, type="l", xlab="time", ylab = "intensity", col="grey")
    lines(est_c$t, est_c$intensity)
    curve(intensity_powerlaw(x,est[1],est[2]), add=TRUE, col="orange")
    curve(intensity_powerlaw(x,a,b), add=TRUE, col="royalblue", lty=2)
  }
}, movie.name = "animation2.gif", interval=0.1)
