library(tidyr)
library(dplyr)
library(ggplot2)

### test functions
LRtest <- function(p,x,n){
  phat <- x/n
  lp <- 2*(dbinom(x, n, phat ,log = TRUE)-
             dbinom(x, n, p ,log = TRUE))
  pchisq(lp, 1, lower.tail = FALSE)
}

Waldtest <- function(p,x,n){
  phat <- x/n
  se <- sqrt(phat*(1-phat)/n)
  2*pnorm(abs(phat-p)/se, lower.tail = FALSE)
}

#Score test is equivalent to prop.test

simfunc <- function(p_true, p_test, n, iter){
    x = rbinom(iter, n, p_true)
    pv_lr = LRtest(p_test, x, n)
    pv_wald = Waldtest(p_test, x, n)
    pv_score = sapply(x, function(x0){prop.test(x0, n, p = p_test, correct = FALSE)$p.value})
    pv_score_c = sapply(x, function(x0){prop.test(x0, n, p = p_test, correct = TRUE)$p.value})
    data.frame(lr=pv_lr, wald=pv_wald, score = pv_score,score_c = pv_score_c, p_true=p_true, p_test=p_test, n=n)
}

sim_settings = expand.grid(p_true = c(0.99, 0.95, 0.9, 0.8, 0.5, 0.2, 0.1, 0.05, 0.01), n=c(10,50,100)) 
dim(sim_settings)
set.seed(1234)
system.time({
  res1 = lapply(seq_len(nrow(sim_settings)), function(i){simfunc(p_true=sim_settings$p_true[i], p_test=sim_settings$p_true[i], n=sim_settings$n[i], 10000)})
})
#  ユーザ システム     経過 
#  21.995    0.057   22.104
# saveRDS(res1, file="result_binom_sim1.rds")
# res1 <- readRDS(file="result_binom_sim1.rds")

res1_df = pivot_longer(bind_rows(res1), lr:score_c,  values_to = "p_value", names_to = "method")

CP_df = group_by(res1_df, method, p_true, n)%>%
  summarise(cp = mean(p_value > 0.05))
CP_df



p_cp = ggplot(CP_df,aes(x=p_true, y=cp, colour=method, linetype=method, shape=method))+
  geom_hline(yintercept =  0.95, colour = "darkgrey")+
  geom_point() + geom_line() +
  scale_y_continuous(n.breaks = 8) +
  labs(y="coverage probability")+
  facet_grid(rows = vars(n), labeller = label_both)+
  scale_color_brewer(palette = "Set2")+
  theme_bw(18) + 
  theme(strip.text.y = element_text(angle=0),strip.background = element_rect(fill = "white"))
print(p_cp)
ggsave("p_cp.png", plot = p_cp, width = 10, height = 10)

p = ggplot(res1_df, aes(colour=method, linetype=method)) + 
  geom_abline(slope = 1, intercept = 0, colour = "lightgrey")+
  stat_ecdf(aes(x = p_value))+
  facet_grid(p_true~n, labeller = label_both)+
  labs(x="nominal", y="actual")+
  scale_color_brewer(palette = "Set2")+
  theme_bw(18) + 
  theme(axis.text.x = element_text(angle = 90),
        strip.text.y = element_text(angle=0), 
        strip.background = element_rect(fill = "white"))
print(p)
ggsave("p_alpha.png", plot = p, width = 10, height = 10)

