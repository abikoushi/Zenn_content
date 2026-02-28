library(deSolve)
library(expm)
#https://cran.r-project.org/web/packages/expm/index.html

modLinear <- function(Time, State, A) {
  return(list(A%*%State))
}

a = 0.2
b = 0.1
A = matrix(c(-a, 0,
             a, -b), 2, 2, byrow = TRUE)
yini  <- c(X = 1, Y = 0)
times <- seq(0, 30, by = 0.5)
out_A   <- ode(yini, times, modLinear, A)

###
# matplot(out_A[,1], out_A[,-1], type = "p", pch=1)
# curve(exp(-a*x), add=TRUE, col=1, lwd=2) #analytic solution
# curve(2*(exp(-b*x)-exp(-a*x)), add=TRUE, col=2, lwd=2)
###

ei_A = eigen(A)
P_A = ei_A$vectors
Pinv_A = solve(ei_A$vectors)

Solm = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

ts = 0:30
sol1 = sapply(ts, Solm, P=P_A, Pinv=Pinv_A, values=ei_A$values, yini=yini)

png("sol1.png", width = 700, height = 700)
matplot(out_A[,1], out_A[,-1], type = "l", xlab = "t", ylab = "U(t)")
matpoints(ts, t(sol1))
dev.off()

###
B = matrix(c(0, 1,
             -1, 0), 2, 2, byrow = TRUE)

out_B   <- ode(yini, times, modLinear, B)
ei_B = eigen(B)
P_B = ei_B$vectors
Pinv_B = solve(P_B)

ts = 0:30
sol2 = sapply(ts, Solm, P=P_B, Pinv = Pinv_B, values=ei_B$values, yini=yini)

print(head(t(sol2)))

png("sol2.png", width = 700, height = 700)
matplot(out_B[,1], out_B[,-1],type = "l", xlab = "t", ylab = "U(t)")
matpoints(ts, t(Re(sol2)))
dev.off()
