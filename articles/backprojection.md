---
title: "Becker et al. (1991) より，感染者数の逆計算を試す"
emoji: "🦠"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: false
---

## 背景

[Becker, N.G., Watson, L.F. and Carlin, J.B. (1991), A method of non-parametric back-projection and its application to aids data. Statist. Med., 10: 1527-1542.](https://doi.org/10.1002/sim.4780101005) の back projection （逆計算）を試してみる．

感染症は潜伏期間があるために発症日がわかっても感染日はわからないことが多い．そこで時刻 $t$ における発症者数 $y_t$ の系列から時刻 $t$ に感染した人の数 $N_t$ の系列を見積もることがしたい．潜伏期間の分布は別の方法によりわかっているものとする．

## モデルとその推定のためのEMアルゴリズム

潜伏期間が $t+1$ である確率を $f_t$ と置く．

$N_t$ の期待値を $\lambda_t$ として、感染がポアソン過程に従っていることにすると、

$$
y_t \sim \operatorname{Poisson}\left(\sum_{i=1}^{t} \lambda_i f_{t-i+1}\right) \tag{1}
$$

となる．

たとえば，

- $t$ 時点に感染して、時刻tに発症する人の数の期待値は $\lambda_t f_1$
- $t-1$ 時点に感染して、時刻tに発症する人の数の期待値は $\lambda_{t-1} f_2$
- $t-2$ 時点に感染して、時刻tに発症する人の数の期待値は $\lambda_{t-2} f_3$
- ...

なので，これを足し合わせたものが $y_t$ の期待値となるということである．

ただし，このままだと対数尤度を計算するときに $\log()$ の中に和が入っている形となってしまい，計算しにくい．そこで生成モデルを次のように改めて考える。

$$
\begin{aligned}
y_t &= \sum_{d=1}^{t} s_{td} \tag{2}\\

\text{where } s_{td} &\sim \operatorname{Poisson}( \lambda_d f_{t-d+1}), \quad d=1, \ldots, t
\end{aligned}
$$

(1)と(2)は等価である．

未観測の $s_{td}$ が得られているとして、対数尤度を書くと，

$$
l(\lambda_t) = \sum_{d=1}^t \{ s_{td} \log(\lambda_t f_{t-d+1}) -\lambda_t f_{t-d+1} \}]
$$

となる．

$\lambda_t$ に関して微分して 0 と置いて解くと，

$$
\hat \lambda_t =\left( \sum_{d=1}^t s_{td} \right) / \left( \sum_{d=1}^t f_{t-d+1}\right)
$$

これが EM アルゴリズムの M ステップである．

E ステップは $s_{td}$ の条件付き期待値を求める．

ポアソン分布の和は和の値で条件つけてやると多項分布になるので，

$$
E[s_{td} | y_t] = y_t \lambda_t f_{t-d+1} / \left(\sum_{d=1}^{t} \lambda_t f_{t-d+1}\right)
$$

となる．

## Rによる数値例

上述した EM アルゴリズムを R で次のように実装した．

```r
lambdaest <- function(count_obs, ft, maxit=1000,tol=0.1){
  len <- length(count_obs)
  lamhat <- rep(mean(count_obs),len)
  den <-rev(cumsum(ft))
  mu <- sapply(1:len, function(t)sum(lamhat[1:t]*f[t:1]))
  ll <- numeric(maxit)
  ll[1] <- sum(dpois(count_obs,mu,log = TRUE))
  for (i in 2:maxit) {
    #E stpp
    yp <- lapply(1:len, function(t)count_obs[t]*lamhat[1:t]*ft[t:1]/sum(lamhat[1:t]*ft[t:1]))
    #M step
    num <- sapply(1:len,function(t)sum(sapply(yp, function(x)x[t]),na.rm = TRUE))
    lamhat <- (num)/(den)
    mu <- sapply(1:len, function(t)sum(lamhat[1:t]*ft[t:1]))
    ll[i] <- sum(dpois(count_obs,mu,log = TRUE))
    if(ll[i]-ll[i-1]<tol){
      break
    }
  }
  return(list(lambda=lamhat,mu=mu,it=i,ll=ll[1:i]))
}
```

シミュレーションでは $t=10,11,12,..20$ では $ \lambda_t = 100$, 他の場合は $\lambda_t = 0$ とした．潜伏期間の分布は形状パラメータ2，スケールパラメータ7のワイブル分布とした．こんなふうだ．

```r
genefun_c <- function(maxrd,maxed,mined,lambda,m,eta){
  np <- rpois(1,(maxed-mined)*lambda)
  ti <- sort(runif(np,mined,maxed))
  si <- rweibull(np,m,eta)
  rd <- ceiling(ti+si)
  id <- ceiling(ti)
  infection <- sapply(1:maxrd,function(t)sum(id==t))
  obs <- sapply(1:maxrd,function(t)sum(rd==t))
  data.frame(obs,infection)
}
```

推定を実行した結果がこちら．まず観測値である発症者数についての当てはまりを確認する．

![](/images/backprojection/fit1.png)
*縦棒が観測された発症者数，丸が推定された平均値*

次に未観測の感染者の推移について確認する．

![](/images/backprojection/est1.png)
*点線が未観測の感染者，丸が推定された平均値*

ちょっとガタガタしてるけどまあまあうまく求まっていることがわかる．

ただし，EMアルゴリズムをどこで止めるかで結果がだいぶかわる．上記のは対数尤度の増加が0.1未満になったら収束と判定した．収束判定を厳しくして閾値を $10^{-8}$ とかにしてみるとこんな感じだ．

![](/images/backprojection/est2.png)
*点線が未観測の感染者，丸が推定された平均値*

尤度は大きくなっていくが $\hat \lambda_t$ はがたがたしてくる．


[Becker et al. (1991)](https://doi.org/10.1002/sim.4780101005) では $\hat \lambda_t$ を事後的にスムージングしている．

これは私見だが，初手では $\hat \lambda_t$ ではなく，次の図のようにその累積和 $\sum_{j=1}^t \hat \lambda_j$ をプロットしてみるのがいいように思う．

![](/images/backprojection/cumulative_intensity.png)

和を取ると大数の法則が効いてくれるので問題が緩和される．これはヒストグラムではビンの幅の取り方に依存して印象が変わるが，経験分布関数ではその煩わしさがないのと類似している．

使用したコード全体は以下：

[https://github.com/abikoushi/Zenn_content/blob/main/R/backprojecttion.R]

おしまい．