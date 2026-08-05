library(dplyr)
library(ggplot2)

simulate_branching <- function(mu,
                               lambda,
                               Tmax) {
  
  N0 <- rpois(1, Tmax*mu)
  if(N0==0L){
    return(NULL)
  }
  t <- sort(runif(N0, 0, Tmax))
  individuals <- data.frame(
    id = 1L:N0,
    parent = 0L,
    birth_time = t,
    generation = 1L
  )
  
  next_id <- N0 + 1L
  i <- 1L
  while (i <= nrow(individuals)) {
    birth <- individuals$birth_time[i]
    t <- birth
    repeat {
      t <- t + rexp(1, rate = lambda)
      if (t > Tmax){
        break
      }
      individuals <- rbind(
        individuals,
        data.frame(
          id = next_id,
          parent = individuals$id[i],
          birth_time = t,
          generation = individuals$generation[i] + 1L
        )
      )
      next_id <- next_id + 1L
    }
    i <- i + 1L
  }
  
  individuals <- individuals[order(individuals$birth_time), ]
  return(individuals)
}

set.seed(805)
tree <- simulate_branching(
  mu = 1,
  lambda = 0.9,
  Tmax = 5
)


segments <- dplyr::filter(tree, parent > 0) |> 
  rowwise() |> 
  mutate(
    x = tree$birth_time[tree$id == parent],
    y = tree$generation[tree$id == parent],
    xend = birth_time,
    yend = generation - 0.05
  )

p <- ggplot(tree) +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.02, "npc")),
    colour="darkgrey") +
  geom_point(aes(birth_time, generation), size=3, shape=1) +
  geom_rug(aes(x=birth_time))+
  scale_y_reverse() +
    labs(x = "time", y="generation") +
    theme_classic()
print(p)
ggsave("tree.png", width = 7, height = 7)

simulate_count <- function(lambda){
  res <- vector("list",100)
  for(i in 1:100){
    tree <- simulate_branching(
      mu = 1,
      lambda = lambda,
      Tmax = 5
    )
    if(length(tree)>0){
      res[[i]] <- data.frame(
        time = tree$birth_time,
        count = seq_along(tree$birth_time),
        group = i
      )
    }
  }
  bind_rows(res)  
}

set.seed(806)
res1 <- simulate_count(0.6)

p <- ggplot(res1, aes(time, count)) +
  geom_step(aes(time, count, group=group), linewidth=0.1) +
  stat_function(fun = function(x)(1/0.6)*expm1(0.6*x), linetype=2, colour="royalblue", linewidth=1)+
  theme_classic()
print(p)
ggsave("count.png", width = 7, height = 7)


totalcounts <- group_by(res1, group) |> 
  summarise(count = last(count)) |> 
  pull(count)

print(mean(totalcounts))
print(var(totalcounts))
