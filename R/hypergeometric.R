library(BiasedUrn)
X <-matrix(c(12,5,6,12), nrow=2)
rs <- rowSums(X)
cs <- colSums(X)
tab <- rbind(tab_drow,tab_all-tab_drow)
#check the marginal
print(all(rowSums(tab)==rs))
print(all(colSums(tab)==cs))


simulated_urn <- function(X, or=1){
  rs <- rowSums(X)
  cs <- colSums(X)
  urn <- rep(factor(c("B","W"), levels = c("B","W")), cs)
  tab_all <- table(urn)
  weight <- rep(c(or,1), cs)
  tab_drow <- table(sample(urn, size = rs[1], prob = weight))
  tab_drow[1]
}


res <- replicate(10000, simulated_urn(X, or=1.5))

simfreq <- table(res)
plot(simfreq/sum(simfreq))
xv <- as.integer(names(simfreq))
points(xv, BiasedUrn::dFNCHypergeo(xv,  cs[1], cs[2], rs[1], 2), type = "b", lty=2)

