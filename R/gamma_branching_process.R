library(dplyr)
library(ggplot2)
library(Rcpp)
Rcpp::sourceCpp("./cpp/gamma_branching_process.cpp")

tree <- simulate_branching_cpp(
  a = 2,
  b = 1,
  Tmax = 10
)

segments <- dplyr::filter(tree, !is.na(parent)) |> 
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


simulate_count <- function(a,b){
  res <- vector("list",5000)
  for(i in 1:5000){
    tree <- simulate_branching_cpp(
      a = a,
      b = b,
      Tmax = 10
    )
    res[[i]] <- data.frame(
      time = c(0, tree$birth_time),
      count = c(0L, seq_along(tree$birth_time)),
      group = i
    )
  }
  bind_rows(res)
}

set.seed(1234)
a <- 2
b <- 0.5
res1 <- simulate_count(a, b)
alpha <- b*(2^(1/a)-1)

p <- ggplot(res1, aes(time, log1p(count))) +
  geom_step(aes(group=group), linewidth=0.01) +
  geom_abline(slope = alpha, colour="orangered", linetype = 2, linewidth = 1.5)+
  theme_classic()
print(p)


totalcounts <- group_by(res1, group) |> 
  summarise(count = last(count)) |> 
  pull(count)

print(mean(totalcounts))
print(expm1(alpha*10))


