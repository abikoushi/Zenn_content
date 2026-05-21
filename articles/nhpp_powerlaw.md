---
title: "べき乗則ポアソン過程の最尤推定"
emoji: "⛓️‍💥"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 確率過程]
published: false
---

## べき乗則ポアソン過程とは

べき乗則ポアソン過程は power-law poisson process の私による直訳で，強度関数（intensity function）が $\lambda(t)= (\beta/\alpha) (t/\alpha)^{\beta-1}$ の非定常ポアソン過程のことである. 信頼性工学などで使われる．

この強度関数はワイブル分布のハザード関数と一致するのでワイブル過程（Weibull process）と呼ばれることもある．ワイブル分布の再生過程と紛らわしいかと思ってこの記事では power-law の方に注目した用語を採用してみた．

パラメータ $\beta$ に応じて「だんだん増える」「だんだん減る」というわかりやすい性質がある．

- $\beta >1$ のとき，事象の起こりやすさ（強度）は時間とともに増加する.
- $\beta=1$ のとき，事象の起こりやすさは時間に依存しない. 点の集合は定常ポアソン過程になる. 点と点の間隔は指数分布に従う.
- $\beta<1$ のとき，事象の起こりやすさは時間とともに減少する.

なんらかのイベント，たとえば故障とか病気とかバグの発見とかが，だんだん増えてるかだんだん減ってるかというのは興味があることが多いだろうと思う．

強度関数 $\lambda(t)$ と累積強度関数

$$
\Lambda(T) = \int_{0}^{T}\lambda\left( u\right) \,du=(t/\alpha)^\beta
$$

をいくつかRを使ってプロットしてみよう．

```r
intensity_powerlaw <- function(t,alpha, beta){
  (beta/alpha)*(t/alpha)^(beta-1)
}

Intensity_powerlaw <- function(x,alpha,beta){
  (x/alpha)^beta
}

#png("intensity.png", width = 700, height=500)
par(mfrow=c(1,2))
curve(intensity_powerlaw(x,1,1), 0, 5, xlab="time", ylab = "intensity", lwd=1.5, ylim=c(0,2))
curve(intensity_powerlaw(x,1,1.5),add=TRUE, col="steelblue",lty=2, lwd=1.5, n=1001)
curve(intensity_powerlaw(x,1,0.8),add=TRUE, col="firebrick",lty=3, lwd=1.5, n=1001)
legend("topright",c("beta=1","beta=1.5","beta=0.8"), 
       lty=1:3, col=c("black","steelblue","firebrick"), lwd=1.5)

curve(Intensity_powerlaw(x,1,1), 0, 5, xlab="time", ylab = "cumulative intensity", lwd=1.5)
curve(Intensity_powerlaw(x,1,1.5),add=TRUE, col="steelblue",lty=2, lwd=1.5)
curve(Intensity_powerlaw(x,1,0.8),add=TRUE, col="firebrick",lty=3, lwd=1.5)
legend("topleft",c("beta=1","beta=1.5","beta=0.8"), 
       lty=1:3, col=c("black","steelblue","firebrick"), lwd=1.5)
#dev.off()
```

![](/images/nhpp_powerlaw/intensity.png)

強度関数は単位時間あたりのイベント発生回数の平均，累積強度関数は累計イベント発生回数の平均に対応する．



## シミュレーションの方法

べき乗則ポアソン過程の擬似乱数を生成する方法を考える．

非定常ポアソン過程にしたがう確率変数 $t_1$ は分布関数 $F(t)=1-\exp(-\Lambda(t))$ を持つ分布に従う．

$z_1, z_2, \ldots, z_n$ を一様乱数として，逆関数法により，$t_1$ は，

$$
t_1=\alpha(- \log (z_1))^{1/\beta}
$$ 

で得られる．つまり最初の一回目のイベントの発生時刻はワイブル分布に従う．

$t_k$ ($k > 1$) の発生時刻については次の条件付き確率で表せる．

$$
 P(T>t_k|T>t_{k-1})= \frac{\exp\left( -( \frac{t_k}{\alpha} )^\beta \right) } {\exp \left(- ( \frac{t_{k-1}}{\alpha})^\beta \right)}
$$

逆関数法により，

$$
t_k=\alpha \left(- \log (z_k) + \left(\frac{t_{k-1}}{\alpha}\right)^\beta  \right)^{1/\beta}．
$$

で得られる．

```r
NHPP_powerlaw = function(Tmax, alpha, beta, maxit){
  mlogz = rexp(1)
  t <- rep(Inf, maxit)
  t[1] = alpha*mlogz^(1/beta[1])
  for(i in 2:maxit){
    mlogz = rexp(1)
    ti = alpha*( mlogz + (t[i-1]/alpha)^beta )^( 1/beta )
    if(ti > Tmax){
      break
    }
    t[i] <- ti
  }
  return(t[1:(i-1)])
}
```

## 最尤推定量

時刻 $t$ における強度が $\lambda(t)$ の非定常ポアソン過程を, 区間 $(0, T]$ で観測した場合の尤度関数は次のようになる.

$$
 L(\alpha,\beta)= \prod_{i=1}^{n} \lambda ( t_{i}) \exp\left(-\Lambda(T)\right)
$$


べき乗則ポアソン過程の場合，パラメータ $\alpha,\beta$ の最尤推定量は閉じた形で求まる．

$$
\begin{aligned}
\log L(\alpha, \beta)&= \left\{ \sum_{i=1}^n \log\lambda ( t_{i})\right\} -\Lambda(T) \\
&= \sum_{i=1}^{n}\left( \log \beta - \log \alpha + (\beta-1)(\log(t_i/\alpha) \right) - (T/\alpha)^\beta \\
&= n\log \beta - n\log \alpha + (\beta-1)\left\{ \sum_{i=1}^{n}\log(t_i/\alpha) \right\}- (T/\alpha)^\beta
\end{aligned}
$$

これを微分して 0 と置いて解く．

$\alpha$ については，

$$
 \frac{\partial}{\partial \alpha} \log L(\alpha, \beta)= -\frac{n}{\alpha} +  \frac{T^\beta}{\alpha^{\beta+1}} = 0,
$$

を解くことにより，

$$
\hat \alpha=\frac{T}{n^{1/\beta}}.
$$

$\beta$ については次の方程式を解く．

$$
 \frac{\partial}{\partial \beta} \log L(\alpha, \beta) = \frac{n}{\beta}+\sum_{i=1}^{n}\log (t_i/\alpha) -(T/\alpha)^\beta \log (T/\alpha)=0.
$$

$\hat \alpha$ を代入して，

$$
\frac{n}{\beta}+\sum_{i=1}^{n}(\log (t_i/T) - (1/\beta)\log n)- (n/\beta) \log n=0,
$$

よって，

$$
\hat \beta=\frac{n}{\sum_{i=1}^{n}\log (T/t_i)}.
$$

以上の議論により，Rでは次のように最尤推定量を実装できる．

```r
estimator <- function(dat, t_max){
  n <-length(dat)
  beta <- n/sum(log(t_max)-log(dat))
  alpha <- t_max/(n^(1/beta))
  c("alpha"=alpha, "beta"=beta)
}
```

## 点推定のシミュレーション

観測期間 $T$ を伸ばしながら最尤推定量の振る舞いを見てみよう．

$\beta=1.5$, $\alpha=1$ のとき：

![](/images/nhpp_powerlaw/animation1.gif)

$\beta=0.8$, $\alpha=1$ のとき：

![](/images/nhpp_powerlaw/animation2.gif)

べき関数にべき関数を当てはめているわけだが観測期間が短いときはけっこうあばれる．観測期間が長くなると収束してくる．

本記事で使用した R のコード全体はこちら：



いったん以上です（区間推定に続くかも）．