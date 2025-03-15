library(animation)
sim_pois_1d <- function(t_n, n, p, seed){
  set.seed(seed)
  counts <- integer(t_n)
  x <- seq(0, 1, length.out=n)
  delta_list <- vector("list", t_n)
  hit <- matrix(0, t_n, n)
  cens = integer(t_n)
  for(i in 1:t_n){
    y <- rbinom(g_n, 1, p)
    counts[i] <- sum(y)
    if(counts[i]>0L){
      hit[i,] <- y
      delta = diff(c(0,x[y==1L])) # geometric
      delta_list[[i]] <- delta
    }
  }
  list(x=x,
       hit=hit,
       delta=delta_list,
       counts=counts)
}

plot_anim <- function(res, rate){
  tab = table(res$counts)
  ran_x = range(res$counts)
  prob_pois = dpois(ran_x[1]:ran_x[2], rate)
  ran_y = c(0, 1)
  bw = 0.1
  bk = seq(0, 1, by=bw)
  animation::saveGIF({
    layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
    for(i in 1:length(res$delta)){
      plot(res$x, res$hit[i,], type = "h", xlab = "x", ylab = "", main=paste("trial:", i))
      tab1 = table(res$counts[1:i])
      plot(tab1/sum(tab1), xlim=ran_x, ylim=ran_y,
           xlab="count", ylab="prob", type = "h", main = "poisson dist.")
      lines(ran_x[1]:ran_x[2], prob_pois, type = "b")
      hist(unlist(res$delta[1:i]), breaks = bk, freq=FALSE,
           main = "exponential dist.", xlab = "diff")
      curve(dexp(x, rate), add=TRUE, lty=2)
    }
  },movie.name = "sim_pois.gif", interval=0.2)
}

n = 100
p = 0.05
rate = g_n*p
res = sim_pois_1d(t_n=100, n=n, p=p, 123)
plot_anim(res, rate)
