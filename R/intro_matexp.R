library(deSolve)
modLinear <- function(Time, State, A) {
  return(list(A%*%State))
}

a = 0.2
b = 0.1
A = matrix(c(-a, 0,
             a, -b), 2, 2, byrow = TRUE)
yini  <- c(X = 1, Y = 0)
times <- seq(0, 30, by = 0.5)
out   <- ode(yini, times, modLinear, A)

matplot(out[,1], out[,-1], type = "p", pch=1)
curve(exp(-a*x), add=TRUE, col=1, lwd=2) #analytic solution
curve(2*(exp(-b*x)-exp(-a*x)), add=TRUE, col=2, lwd=2)

ei_A = eigen(A)
ei_A$vectors
P = ei_A$vectors
Pinv = solve(ei_A$vectors)

SolA = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

ts = 0:30
sol0 = sapply(ts, SolA, P=P, Pinv=Pinv, values=ei_A$values, yini=yini)
matplot(out[,1], out[,-1], type = "l")
matpoints(ts, t(sol0))
