library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(MASS)

make_prob <- function(alpha, eta){
  xi <- plogis(outer(-eta, alpha, "+"))
  for(j in 2:length(alpha)){
    xi[,j] <- xi[,j]*(1-rowSums(xi[,1:(j-1), drop=FALSE]))  
  }
  cbind(xi, 1-rowSums(xi))  
}

make_prob_l <- function(alpha, eta){
  xi <- plogis(outer(-eta, alpha, "+"))
  for(j in 2:length(alpha)){
    xi[,j] <- xi[,j]-xi[,j-1]
    }
  cbind(xi, 1-rowSums(xi))  
}

make_prob_p <- function(alpha, eta){
  xi <- pnorm(outer(-eta, alpha, "+"))
  for(j in 2:length(alpha)){
    xi[,j] <- xi[,j]-xi[,j-1]
  }
  cbind(xi, 1-rowSums(xi))  
}


PHreg <-function(Y, X, prec = 0, method = "BFGS"){
  K <- ncol(X)
  J <- ncol(Y)
  M <- t(apply(Y,1, function(x)rev(cumsum(rev(x)))))
  LL <- function(par, Y, M, J, X, K, prec){
    alpha <- par[1:(J-1)]
    beta <- par[J:(J+K-1)]
    out <- 0
    Xbeta <- -X%*%beta 
    for(j in 1:(J-1)){
      out <- out +
        sum(Y[,j]*(alpha[j]+Xbeta) -
              M[,j]*log1p(exp(alpha[j]+Xbeta)))
    }
    out <- out - 0.5*sum(beta^2*prec)
    return(out)
  }
  
  dLL <- function(par, Y, M, J, X, K, prec){
    alpha <- par[1:(J-1)]
    beta <- par[J:(J+K-1)]
    dLLalpha <- function(alpha,beta,Y,M,J,X){
      Xbeta <- -X%*%beta
      out <- numeric(J-1)
      for(j in 1:(J-1)){
        out[j] <- sum(Y[,j] - M[,j]*plogis(alpha[j]+Xbeta))
      }
      return(out)
    }
    
    dLLbeta <- function(alpha, beta, Y, M, J, X, K, prec){
      Xbeta <- -drop(X%*%beta)
      out <- numeric(K)
      for(j in 1:(J-1)){
        out <- out + colSums(-Y[,j]*X + M[,j]*X*plogis(alpha[j]+Xbeta))
      }
      out <- out - beta*prec
      return(out)
    }
    
    return(
      c(dLLalpha(alpha, beta, Y, M, J, X),
        dLLbeta(alpha, beta, Y, M, J, X, K, prec))
    )
  }
  ini <-numeric(J+K-1)
  xnames <- colnames(X)
  if(is.null(xnames)){
    xnames <- paste0("X", seq_len(ncol(X)))
  }
  names(ini) <- c(paste0("alpha", 1:(J-1)), xnames)
  opt <- optim(ini, LL, dLL, method = method,
               control = list(fnscale=-1),
               J=J, Y=Y, M=M, X=X, K=K, prec=prec, 
               hessian = TRUE)
  return(opt)
}



fit_l <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, 
              method = "logistic", Hess = TRUE)

CI_l <- confint(fit_l)
colnames(CI_l) <- c("lower", "upper")
CI_l <- rownames_to_column(
    data.frame(CI_l, est=fit_l$coefficients, method="logistic")
  )

fit_p <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing, 
              method = "probit", Hess = TRUE)
CI_p <- confint(fit_p)
colnames(CI_p) <- c("lower", "upper")
CI_p <- rownames_to_column(
    data.frame(CI_p, est=fit_p$coefficients, method="probit")
  )

housing_wide <- pivot_wider(housing, names_from = Sat, values_from = Freq)
Y <- as.matrix(housing_wide[,4:6])
X <- model.matrix(~ Infl+Type+Cont, housing_wide)
X <- X[,-1]

res <- PHreg(Y,X)
res$value

####
#goodeness-of-fit
housing_wide_prob <- group_by(housing,  Infl, Type, Cont) %>% 
  mutate(prob = Freq/sum(Freq)) %>%
  ungroup() %>% 
  dplyr::select(-Freq) %>% 
  pivot_wider(names_from = Sat, values_from = prob)

Prob_ph <- make_prob(alpha = unname(res$par[1:2]), eta = c(X%*%res$par[3:8]))
Prob_l <- make_prob_l(alpha = unname(fit_l$zeta), eta= c(X%*%fit_l$coefficients))
Prob_p <- make_prob_p(alpha = unname(fit_p$zeta), eta= c(X%*%fit_p$coefficients))

dffit <- data.frame(Prob_obs) %>% 
  pivot_longer(1:3, names_to = "Sat", values_to = "obs") %>% 
  mutate(PH = c(t(Prob_ph)),
         logistic = c(t(Prob_l)),
         probit = c(t(Prob_p)),
         Sat = factor(Sat, levels = c("Low", "Medium", "High"))) %>% 
  pivot_longer(PH:probit, names_to = "method", values_to = "fit") 

ggplot(dffit, aes(x=obs, y=fit, colour = method, shape=method))+
  geom_abline(slope = 1, linetype=2)+
  geom_point(alpha=0.75)+
  facet_wrap(~Sat)+
  theme_bw()

ggsave("fit.png", width = 5, height = 3)
# sum(as.matrix(housing_wide[,4:6])*log(Prob_ph))
# sum(as.matrix(housing_wide[,4:6])*log(Prob_p))
# sum(as.matrix(housing_wide[,4:6])*log(Prob_l))

print(res$value)
print(logLik(fit_l))
print(logLik(fit_p))

#####
#CI

se <- sqrt(-diag(solve(res$hessian)))
q <- qnorm(0.975)
CI_ph <- data.frame(lower=res$par - q*se,
                    upper=res$par + q*se,
                    est=res$par, method="PH")
CI_ph <- rownames_to_column(CI_ph[3:8,])

CI <- bind_rows(CI_p, CI_l, CI_ph)

ggplot(CI, aes(x= rowname, y=est, colour=method, shape=method))+
  geom_pointrange(aes(ymin=lower,ymax=upper), position = position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype=2)+
  labs(x="variable", y="estimates")+
  coord_flip()+
  theme_bw(12)

ggsave(filename = "coef.png", width = 6, height = 4)

