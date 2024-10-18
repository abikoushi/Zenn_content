library(BiasedUrn)
library(exact2x2)
library(animation)

#MCMCpack::dnoncenhypergeom(x = NA, cs[1],cs[2],rs[1], 1)
saveGIF({
  for(i in c(0:17,17:0)){
    rs <- c(18,17)
    cs <- c(17,18)
    X <-matrix(c(i,cs[1]-i,rs[1]-i,cs[2]-(rs[1]-i)), nrow=2)
    pfun_f <- function(phi){
      sapply(phi, function(phi){exact2x2(X, or=exp(phi))$p.value})
    }
    lp_f <-  function(phi){
      sapply(phi, function(phi){BiasedUrn::dFNCHypergeo(X[1,1], 17,cs[2],18, exp(phi))})
    }
    #range (prior distribution's support)
    mx <- 4
    mn <- -4
    int <- integrate(lp_f, mn, mx) #normalize constant
    postfun_F <- Vectorize(function(u){
      lower <- integrate(lp_f, mn, u)$value/int$value
      upper <- integrate(lp_f, u, mx)$value/int$value
      2*pmin(lower, upper)
    })
    curve(pfun_f(x), from = mn, to =mx, lty=2, ylab = "p-value", n=501, ylim=c(0,1))
    curve(postfun_F, add=TRUE, n=501, col="royalblue")
  }
},  movie.name = "fisher.gif", interval=0.2)
# legend("topleft", c("fisher's test","Bayesian"),
#        lty=2:1, col = c("black","royalblue"), box.lwd = 0)    
