iter_nextstep <- function(n, K){
  I <- c(1, 0)
  I_hist <- matrix(0, nrow = 2, ncol = n+1)
  I_hist[, 1] <- I
  for(i in 1:n){
    I <- K%*%I
    I_hist[,i+1] <- I  
  }
  return(t(I_hist))  
}

K <- matrix(c(1.4, 0.4,
              0.6, 0.3), byrow = TRUE, nrow = 2)
n <- 10
I_hist <- iter_nextstep(n,K)

png("original1.png")
matplot(0:n, I_hist, type = "b", pch = c("C","A"),
        xlab = "t", ylab = expression(I[t]))
dev.off()

ei_K <- eigen(K)
print(ei_K$values)

norm_v <- ei_K$vectors[,1]/sum(ei_K$vectors[,1])
print(norm_v)

sumI <- rowSums(I_hist)
normI <- sweep(I_hist, 1, sumI,"/")

png("ratio11.png")
plot(sumI[-1]/sumI[-length(sumI)], type = "b",
     ylab = expression(N[t+1]/N[t]), xlab = "t")
abline(h=ei_K$values[1],lty=3)
dev.off()

png("ratio21.png")
matplot(0:n, normI, type = "b",ylab = expression(I[t]/N[t]), xlab = "t",
        pch = c("C","A"))
abline(h = norm_v, lty=3)
dev.off()

###

K2 <- matrix(c(0, 2,
              3, 0), byrow = TRUE, nrow = 2)

ei_K2 <- eigen(K2)
print(ei_K2$values)

I_hist <- iter_nextstep(n,K2)

png("original2.png")
matplot(0:n, I_hist, type = "b",pch = c("C","A"))
dev.off()

sumI <- rowSums(I_hist)
normI <- sweep(I_hist,1,sumI,"/")

png("ratio12.png")
plot(sumI[-1]/sumI[-length(sumI)], type = "b",
     ylab = expression(N[t+1]/N[t]), xlab = "t")
dev.off()

png("ratio22.png")
matplot(0:n, normI, type = "b",
        ylab = expression(I[t]/N[t]), xlab = "t", pch = c("C","A"))
dev.off()
