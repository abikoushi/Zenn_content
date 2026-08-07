library(dplyr)
library(ggplot2)
?dnbinom
###
branching_process <- function(R0, k, G){
  gen <- vector("list", G)
  p <- k/(k+R0)
  X <- rnbinom(1, k, p) #第1世代
  if(X>0){
    id <- 1L:X
    parent <- 0L
    gen[[1]] <- data.frame(id = id,
                           parent = parent,
                           g = 1L)
    parent <- id
    next_id <- X
    
    for(i in 2L:G){  #第2世代以降
      X <- rnbinom(length(parent), k, p)
      sumX <- sum(X)
      if(sumX == 0){ break } #感染者が0になったら終了
        parent <- rep(parent, X)
        id <- (next_id+1L):(next_id+sumX)
        gen[[i]] <- data.frame(id = id,
                               parent = parent,
                               g = i)
        parent <- id
        next_id <- next_id+sumX
    }    
  }
  dplyr::bind_rows(gen)
}


prob_totalsize <- function(y,R0,k){
  if(y==1){
   1/(1+R0/k)^k
  }else{
    prod(0:(y-2L)/k+y)*((k/(R0+k))^(k*y))*(R0*k/(R0+k))^(y-1)/factorial(y) 
  }
}


rep_sim_branch <- function(iter,R0,k,G){
  res <- vector("list",iter)
  total_size <- integer(iter)
  for( i in 1:iter){
    tree <- branching_process(R0, k, G)
    if(length(tree)>0){
      tree$group <- i
      res[[i]] <- tree
      total_size[i] <- nrow(tree)
    }else{
      res[[i]] <- data.frame(id=NA_integer_,
                             parent=NA_integer_,
                             g=NA_integer_,
                             group=i)
      total_size[i] <- 0L
    }
  }
  res <- bind_rows(res)
  list(tree=res, total_size=total_size)
}

set.seed(12345)
k <- 5
R0 <- 0.6
G <- 20
# tree <- branching_process(R0, k, G)
res <- rep_sim_branch(10000,R0,k,G)

tab <- table(res$total_size)
p <- sapply(1:32, prob_totalsize, R0=R0,k=k)
plot(tab/sum(tab), ylim = c(0, 0.6))
points(0:31, p)

max(res$tree$g, na.rm = TRUE)

ggplot(res)+
  geom_point(aes(x=id, y=g), shape=1, size=3)+
  geom_segment(aes(yend=g, xend=id, x = parent, y=g-0.95), 
               arrow = arrow(length = unit(0.2,"cm")))+
  scale_y_reverse()+
  labs(x = "node id", y = "generation")+
  theme_classic()+
  facet_wrap(~group)
