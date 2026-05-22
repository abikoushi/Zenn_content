---
title: "ガンマ分布の形状パラメータの検定（あるいは局外母数があるときのスコア検定の例題）"
emoji: "📕"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---

## はじめに

ガンマ分布の密度関数は次式で与えられる．

$$
f(x;\alpha,\theta)= \frac{1}{\Gamma(\alpha)\theta^\alpha} x^{\alpha-1}e^{-x/\theta}, \qquad x>0.
$$

$\alpha$ を形状（shape）パラメータ， $\theta$ を尺度（scale）パラメータと呼ぶ．

ガンマ分布は形状パラメータによって大きく性質が変わる. 形状パラメータが 1 より大きいときはモードの周りに集中した分布になり，1 以下のときは 0 で発散する裾の重い分布になる.

仮にお金の分布だとすると形状パラメータが 1 以下のときは貧富の差が大きい社会と言えそうだ．

一方で尺度パラメータが $k$ 倍になることはガンマ分布に従う確率変数を $k$ 倍することと等しい．尺度パラメータは文字通りスケールの取り方に依存して変わりうる．

そこで尺度パラメータを局外パラメータ（nuisance parameter）として形状パラメータの検定を行うことを考えてみる．

すなわち，$X_1,\dots,X_n$ を独立に形状パラメータ $\alpha>0$，尺度パラメータ $\theta>0$ のガンマ分布に従う確率変数とし，このときに，

- 興味のあるパラメータ：形状パラメータ ($\alpha$)
- 局外パラメータ：尺度パラメータ ($\theta$)

として，帰無仮説 $\alpha=\alpha_0$ を検定するスコア検定を導く．

## 検定統計量の導出

### 最尤推定量

対数尤度は

$$
L(\alpha,\theta)=(\alpha-1)\sum_{i=1}^n \log X_i-\frac1\theta\sum_{i=1}^n X_i-n\alpha\log\theta-n\log\Gamma(\alpha).
$$

形状パラメータ $\alpha$ に関するスコア関数は

$$
U_\alpha(\alpha,\theta)=
 \frac{\partial L}{\partial \alpha}L(\alpha,\theta)=
 \left\{\sum_{i=1}^n \log X_i \right\}- n\log\theta- n\psi(\alpha),
$$

ここで $\psi(\alpha)=\Gamma'(\alpha)/\Gamma(\alpha)$ はディガンマ関数．

尺度パラメータ $\theta$ に関するスコア関数は，

$$
U_\theta(\alpha,\hat\theta)= \frac{\partial}{\partial \theta} L(\alpha,\hat\theta)= 
 \frac{1}{\theta^2}\sum_{i=1}^n X_i -\frac{n\alpha}{\theta}.
$$


帰無仮説 $\alpha=\alpha_0$ の下で $\theta$ の最尤推定量を求めると，$U_\theta=0$ を解くことにより，

$$
\hat\theta_0 =
\frac{\bar X}{\alpha_0},
$$

ただし

$$
\bar X=\frac1n\sum_{i=1}^n X_i.
$$


$\hat\theta_0$ を $U_\alpha$ に代入すると，
$$
U_\alpha(\alpha_0,\hat\theta_0)=n\left(\overline{\log X}-\log \bar X+\log \alpha_0-\psi(\alpha_0)\right)
$$

ただし

$$
\overline{\log X} = \frac1n\sum_{i=1}^n \log X_i.
$$

### 条件付きスコア検定

スコア検定は，スコア関数がフィッシャー情報行列の逆行列を分散共分散行列とした多変量正規分布に漸近的に従うことを利用したものであった．フィッシャー情報行列の成分は次のようになる．

$$
I_{\alpha\alpha}=n\psi'(\alpha),
$$

$$
I_{\alpha\theta}=\frac{n}{\theta},
$$

$$
I_{\theta\theta}=\frac{n\alpha}{\theta^2},
$$

ここで $\psi'(\alpha)$ はトリガンマ関数である．

$\theta$ のスコアを所与としたときの条件付き分散を求めてみる．

$$
I_{\mathrm{eff}} = I_{\alpha\alpha} - I_{\alpha\theta} I_{\theta\theta}^{-1}
I_{\theta\alpha}
=
n\left(\psi'(\alpha)-\frac1\alpha\right).
$$

これを有効情報量（efficient information）と呼ぶことがある．$\theta$ のスコアを所与としたとき次の検定統計量 $T$ は漸近的に標準正規分布に従う．

$$
T=\frac{U_\alpha(\alpha_0,\hat\theta_0)}{\sqrt{n\left(\psi'(\alpha_0)-\frac1{\alpha_0}\right)}}.
$$

あるいは $T^2$ が自由度1のカイ二乗分布に従うことを利用してもよい．


## R によるシミュレーション

両側検定の検定統計量と p 値は，次のように書ける．

```r
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
```

まず帰無仮説の $\alpha_0$ を真値の $\alpha$ と同じ値に設定し，サンプルサイズ $n$ と $\alpha$ を動かして第一種の過誤の確率をシミュレーションしてみる．次のような関数を書いた．

```r
pv_simfun <- function(n, shape0, shape, scale){
  pv <- numeric(10000)
  for(i in 1:10000){
    x <- rgamma(n, shape = shape, scale = scale)
    pv[i] <- gamma_shape_score_test(x, shape0)$p.value  
  }
  return(pv)
}
```

![](/images/gam_shp_score/alpha.png)

おおむね名目上の水準が保たれていることがわかる．

次に $\alpha_0 = 1$ として検出力のシミュレーション．

![](/images/gam_shp_score/beta.png)

とりあえず導出したスコア検定が検定として使えそうなことはわかった．

