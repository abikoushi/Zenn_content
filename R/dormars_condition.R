library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)
#data are available in:
#https://www.esri.cao.go.jp/jp/sna/menu.html
#https://japanyieldcurve.com/data/

# ファイルのパスは書き換えてください
gdp <- read_csv("~/data/dormar/ritu-mfy2622.csv", skip = 6)
gdp <- mutate(gdp, year = as.integer(gsub("/4-3\\.","",`Fiscal Year`))) |> 
  rename(g = ...2)

# ファイルのパスは書き換えてください
yields <- read_csv("~/data/dormar/jgb-yields-history.csv")
yields_sum <-rename(yields, r = `10Y`) |> 
  dplyr::filter(!is.na(r)) |> 
  mutate(year = as.integer(substr(as.character(date),1,4))) |> 
  group_by(year) |> 
  summarise(mu = mean(r), sigma=sd(r))


ggplot(yields_sum,aes(x=year))+
  geom_line(aes(y=mu, linetype = "r"), colour="orange")+
  geom_ribbon(aes(ymax=mu+sigma, ymin=mu-sigma), fill="orange", alpha=0.1)+
  geom_line(data = gdp, aes(y=g, linetype = "g"),colour="royalblue")+
  labs(linetype="", y="日本の10年国債の金利 r と名目GDP成長率 g（％）", x="年")+
  theme_classic(base_size = 12, base_family = "Osaka")

ggsave("gr_japan.png", width = 8, height = 8)

####

debt_seq <- function(t, r, g, d, b1){
  b <- numeric(t)
  b[1] <- b1
  A <- (1+r)/(1+g)
  for(i in 1:(t-1)){
    b[i+1] <- A*b[i] + d      
  }
  data.frame(t=1:t, b=b)
}

debt_n <-function(n, r, g, d, b1){
  A <- (1+r)/(1+g)
  A^(n-1) * b1 + d*((1-A^(n-1))/(1-A))
}

res1 <- debt_seq(500, r=0.01, g=0.02, d=0.02, b1=1)
res1_2 <- sapply(1:500, function(x)debt_n(x, r=0.01, g=0.02, d=0.02, b1=1))

plot(res1, type = "l")
lines(res1_2, col="red", lty=2)

res2 <- debt_seq(500, r=0.03, g=0.015, d=-0.01, b1=1)
res2_2 <- sapply(1:500, function(x)debt_n(x, r=0.03, g=0.015, d=-0.01, b1=1))

plot(res2, type = "l")
lines(res2_2, col="red", lty=2)

p1 <- ggplot(res1,aes(x=t,y=b))+
  geom_line()+
  labs(title="r=0.01, g=0.02, d=0.02")+
  theme_classic(base_size = 12, base_family = "Osaka")

p2 <- ggplot(res2,aes(x=t,y=b))+
  geom_line()+
  labs(title="r=0.03, g=0.015, d=-0.01")+
  theme_classic(base_size = 12, base_family = "Osaka")

p1/p2

ggsave("Fig31.png", width = 8, height = 8)
