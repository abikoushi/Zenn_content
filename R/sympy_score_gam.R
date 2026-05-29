gamma_shape_score_test <- function(x, alpha0) {
  n <- length(x)
  xbar <- mean(x)
  logxbar <- mean(log(x))
  
  U <- n * (logxbar - log(xbar) + log(alpha0) - digamma(alpha0))
  Ieff <- n * (trigamma(alpha0) - 1 / alpha0)
  stat <- (U^2) / Ieff
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  
  list(
    statistic = stat,
    p.value = pval,
    null.value = alpha0
  )
}

cols <- c("orange","grey20","royalblue")
png("density.png")
curve(dgamma(x,0.5),xlim = c(0,6), col=cols[1], lty=1, lwd=1.5, ylab = "density")
curve(dgamma(x,1),add=TRUE,col=cols[2], lty=2, lwd=1.5)
curve(dgamma(x,2),add=TRUE,col=cols[3], lty=3, lwd=2)
legend("topright", paste("alpha =",c("0.5","1","2")), lty=1:3, col=cols, lwd=2)
dev.off()

####
library(reticulate)
#virtualenv_install("r-reticulate", "sympy")
reticulate::py_require("sympy")
sympy <- import("sympy")
stats <- import("sympy.stats")


alpha <- sympy$Symbol("alpha", positive = TRUE)
theta <- sympy$Symbol("theta", positive = TRUE)
X <- stats$Gamma("X", alpha, theta)

logpdf <- sympy$simplify(sympy$log(stats$density(X)(X)))

# score
S_alpha <- sympy$simplify(sympy$diff(logpdf, alpha))
S_theta <- sympy$simplify(sympy$diff(logpdf, theta))

# Fisher information
I_11 <- sympy$simplify(-stats$Expectation(sympy$diff(S_alpha, alpha)))
I_12 <- sympy$simplify(-stats$Expectation(sympy$diff(S_alpha, theta)))
I_21 <- sympy$simplify(-stats$Expectation(sympy$diff(S_theta, alpha)))
I_22 <- sympy$simplify(-stats$Expectation(sympy$diff(S_theta, theta)))

I_eff <- sympy$simplify((I_11 - I_12 * I_21 / I_22))

print(I_eff)

# データ
set.seed(1234)
x_data <- rgamma(5, shape = 2)

#帰無仮説のalpha
alpha0 <- 1

# スコアの和
S_alpha0 <- sympy$Integer(0)
S_theta0 <- sympy$Integer(0)
I_eff0  <- sympy$Integer(0)
for(val in x_data) {
  S_alpha0 <- S_alpha0 + S_alpha$subs(dict(alpha = alpha0, X = val))
  S_theta0 <- S_theta0 + S_theta$subs(dict(alpha = alpha0, X = val))
  I_eff0 <- I_eff0 + I_eff$subs(dict(alpha=alpha0))
}


sol_theta <- sympy$solve(sympy$Eq(S_theta0, 0), theta)

U_eff0 <- sympy$simplify(
  S_alpha0$subs(dict(theta=sol_theta[[1]]))^2/I_eff0$subs(dict(theta=sol_theta[[1]]))
)


builtins <- import_builtins()
U_eff0_val <- py_to_r(builtins$float(U_eff0$evalf()))

res_r <- gamma_shape_score_test(x = x_data, 1)

print(U_eff0_val)
print(res_r$statistic)
print(all.equal(res_r$statistic, U_eff0_val))


