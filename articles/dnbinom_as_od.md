---
title: "ポアソン分布より分散が大きい分布としての負の二項分布"
emoji: "🥚"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計]
published: true
---

## 基礎

ポアソン分布はカウントデータ（非負の整数のデータ）に対する基本的な分布で、平均と分散が等しい性質があります。

一方、現実のデータ分析では平均より分散が大きいケースがよくあり、そのような状況を過分散（overdispersion）と呼んだりします。

過分散に対応するため、ポアソン分布の代わりに負の二項分布が使われることがあります。
このことはパラメータがガンマ分布に従って変化するようなポアソン分布が負の二項分布になることを知ると理解しやすいです。

ポアソン分布のパラメータを $r \mu$ として、$\mu$ は固定で $r$ がガンマ分布に従って変化する場合、どのような分布になるかを考えます。
すなわち、次のようなデータ生成過程を考えます。

$$
\begin{align*}
y | r &\sim \mathrm{Poisson}(y|r \mu) \\
r &\sim \mathrm{Gamma}(r|\alpha, \alpha)
\end{align*}
$$

ここではガンマ分布の2つのパラメータをどちらも同じ $\alpha$ として平均を1に固定しています。

$r$ を積分消去します。

$$
\begin{align*}
p(y) &= \int p(y|r)p(r) dr\\
&= \int^{\infty}_{0} \mathrm{Poisson}(y|r \mu) \cdot \mathrm{Gamma}(r|\alpha, \alpha) dr\\
&=\int^{\infty}_{0} \frac{(r \mu)^y}{y!}e^{-r \mu} \times \frac{\alpha^{\alpha}}{\Gamma(\alpha)}r^{\alpha-1}e^{-\alpha r}\,dr
\end{align*}
$$

$r$ に依存しない因子を積分記号の外に出して計算します。

$$
\begin{align*}
p(y)&=\frac{\mu^y \alpha^\alpha}{y!\Gamma(\alpha)}\left(\frac{1}{\mu+\alpha}\right)^{y+\alpha}\int^{\infty}_{0} r^{y+\alpha}e^{-r}\,dr\\
&=\frac{\Gamma(y+\alpha)}{y!\Gamma(\alpha)}\left(\frac{1}{\mu+\alpha}\right)^{y+\alpha}\mu^y\alpha^\alpha\\
&=\frac{\Gamma(y+\alpha)}{y!\Gamma(\alpha)}\left(\frac{\alpha}{\alpha+\mu}\right)^\alpha\left(\frac{\mu}{\alpha+\mu}\right)^y
\end{align*}
$$

これは負の二項分布 $\mathrm{NegBin}(\alpha, \mu)$ です。


統計ソフト R の負の二項分布関連の関数（`dnbinom`、`pnbinom`、`qnbinom`、`rnbinom`) では $\alpha$ は `size` という引数、 $\mu$ は `mu` という引数で指定します。

念のため、ガンマ分布のレート（スケールの逆数）パラメータが自由に動くときも考えておきます。

$$
\begin{align}
y | r &\sim \mathrm{Poisson}(y|r \mu) \\
r &\sim \mathrm{Gamma}(r|\alpha, \beta)
\end{align}
$$

ガンマ分布のレートパラメータを $\beta$ 倍することはガンマ分布に従う確率変数を $1/\beta$ 倍することと同値なので、(1)-(2) のモデルで表される $y$ は次の (3)-(4) で表される $y$ と同値です。

$$
\begin{align}
y | r &\sim \mathrm{Poisson}\left(y| r  \cdot (\alpha/\beta) \cdot \mu\right) \\
r &\sim \mathrm{Gamma}(r|\alpha, \alpha)
\end{align}
$$

$\mu' = (\alpha/\beta) \cdot \mu$ とまとめると、上とまったく同じ議論により、やはり負の二項分布が得られます。

$$
p(y)=\frac{\Gamma(y+\alpha)}{y!\Gamma(\alpha)}\left(\frac{\alpha}{\alpha+\mu'}\right)^\alpha\left(\frac{\mu'}{\alpha+\mu'}\right)^y
$$

## R によるケーススタディ

[統計局ホームページ)/第六十三回日本統計年鑑 平成26年−第2章 人口・世帯 ](https://www.stat.go.jp/data/nenkan/back63/02.htm) の「2 - 5 都道府県別人口集中地区人口，面積及び人口密度（エクセル：34KB）」 を使い、人口の分布を調べます。

```r
library(dplyr)
library(ggplot2)
library(readxl)
library(ggrepel)
dat <-read_xls("~/Downloads/y0205000.xls",skip=8) %>% 
  slice(-c(1:8)) %>% 
  select('...1', '人口...3','面積...4') %>% 
  setNames(c("pref","pop","area")) %>% 
  dplyr::filter(!is.na(pop)) %>% 
  mutate(pop=as.integer(pop),area=as.numeric(area)) %>% 
  mutate(pref=gsub("^[0-9][0-9]　", "", pref))
```

まずポアソン分布を用いたモデルを考え、「面積あたりの人口」をパラメータとします。

すなわち、次のような確率分布を考えます。

$$
y \sim \mathrm{Poisson}(y|\tau \lambda)
$$

ここで $y$ を人口、$\tau$ を面積としました。

最尤推定量は次のように得られます（証明略）。

```r
lambahat <- sum(dat$pop)/sum(dat$area)
```

最尤推定されたパラメータを持つ分布から1万回乱数の抽出をおこない、分布の95%区間をプロットしてみます。

```r
set.seed(1234)
sim1 <-matrix(rpois(10000*47,dat$area*lambahat),47,10000)
int1 <-as.data.frame(t(apply(sim1,1,quantile,prob=c(0.025,0.975)))) %>% 
  setNames(c("lower","upper")) %>% 
  mutate(area=dat$area)

ggplot(dat,aes(x=area))+
  geom_point(aes(y=pop))+
  geom_ribbon(data=int1,aes(ymin=lower,ymax=upper),alpha=0.3)+
  theme_bw(16)
```

![](/images/dnbinom_as_od/pop_pois.png)


縦軸が人工、横軸が面積です。

平均はだいたい捉えていますが、データの散らばり具合をカバーしきれていないことがわかります。

次に負の二項分布を用いて、

$$
y \sim \mathrm{NegBin}(y|\alpha,\tau \mu)
$$

なるモデルを考えます。

最尤推定は最適化関数の `optim` に丸投げします。

```r
ll <- function(par,y,tau){
  sum(dnbinom(y,size = exp(par[1]), mu = exp(par[2]+tau), log = TRUE))
}

opt <-optim(c(0,0),ll,y=dat$pop,tau=log(dat$area),control = list(fnscale=-1))
```
こちらも最尤推定されたパラメータを持つ分布から1万回乱数の抽出をおこない、分布の95%区間をプロットしてみます。

前半の考察を確かめるため、ガンマ分布にしたがってパラメータが変化するポアソン分布を愚直に書いていdます。

```r
r <- rgamma(10000*47,exp(opt$par[1]),exp(opt$par[1]))
sim2 <-matrix(rpois(10000*47,dat$area*r*exp(opt$par[2])),47,10000)
int2 <-as.data.frame(t(apply(sim2,1,quantile,prob=c(0.025,0.975)))) %>% 
  setNames(c("lower","upper")) %>% 
  mutate(area=dat$area)


subdat <- dat[dat$pop > int2$upper,]

ggplot(dat,aes(x=area))+
  geom_point(aes(y=pop))+
  geom_ribbon(data=int2,aes(ymin=lower,ymax=upper),alpha=0.3)+
  geom_text_repel(data = subdat , aes(y=pop, label = pref), family="Osaka")+
  theme_bw(16)

```

![](/images/dnbinom_as_od/pop_negbin.png)

データの散らばり具合に近い分布が推定されていることがわかります。