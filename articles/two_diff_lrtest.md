---
title: "2つの尤度比検定：正規分布による例題"
emoji: "🍒"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: false
---

## はじめに

ネイマン・ピアソンの補題により最強力検定になる尤度比検定は帰無仮説と対立仮説がともに一点で与えられるときの尤度比検定であり，最尤推定量を用いる尤度比検定とは異なるので注意が必要かもしれない．ここでは簡単な例題をやってみる．

例題の設定は [竹村彰通『現代数理統計学』（学術図書出版社）](https://www.gakujutsu.co.jp/product/978-4-7806-0860-1/) 第8章を参考にした．

## 手計算してみる

### ネイマン・ピアソンの補題により最強力検定になる尤度比検定

平均 $\mu$ が未知，分散1が既知の正規分布に従うランダム標本 $X_i$ ($i=1, \ldots, n$) を考える．

$$
X_i \sim \mathcal{N}(\mu,1).
$$

（実は，平均が未知だが分散がわかっている状況というのはあまりない．だから人工的な例で恐縮だが，この例で後述の最尤推定量を用いる尤度比検定との比較はできる．）

帰無仮説 $H_0 : \mu = \mu_0$，対立仮説 $H_1 : \mu = \mu_1$ の統計的仮説検定を考える．

尤度比 $L(\mu_1, \mu_0)$ は次のように書ける．

$$
\begin{aligned}
L(\mu_1, \mu_0) &= \exp\left( \sum_{i=1}^{n}\frac{1}{2}(X_i-\mu_1)^2 - \sum_{i=1}^{n}\frac{1}{2}(X_i-\mu_0)^2 \right) \\
&= \exp\left( (\mu_1 - \mu_0)\left(\sum_{i=1}^{n}X_i\right)- \frac{n}{2}(\mu_1^2 - \mu_0^2) \right).
\end{aligned}
$$

尤度比検定ではこの尤度比 $L(\mu_1, \mu_0)$ がある定数 $C$ より大きいときに帰無仮説を棄却する．$L(\mu_1, \mu_0) > C$ となることは，

$$
\frac{1}{n}\sum_{i=1}^{n}X_i> \frac{\frac{1}{n}\log C + \frac{1}{2}(\mu_1^2 - \mu_0^2)}{\mu_1 - \mu_0}
$$

と同値である．左辺を $\bar{X}$，右辺を $C'$ と置く．いま，帰無仮説の下で評価した確率 $P(\bar{X}>C')$ を所与の有意水準 $\alpha$ と一致させたい．$\sqrt{n}(\bar{X} - \mu_0)$ が帰無仮説の下で標準正規分布に従うため，

$$
P(\bar{X}>C') = P(\sqrt{n}(\bar{X} - \mu_0)>\sqrt{n}(C'-\mu_0))
$$

より，$z_\alpha = \sqrt{n}(C'-\mu_0)$ とすればよい．ここで $z_\alpha$ は標準正規分布の上側 $\alpha$ 点とした．すなわち，

$$
C' = \mu_0 + \frac{1}{\sqrt{n}}z_\alpha .
$$

また，P 値は帰無仮説を棄却できる最小の有意水準であるので，標準正規分布の分布関数 $\Phi(x)$ を用いて，次で得られる．

$$
\text{P-value} = 1-\Phi(\sqrt{n}(\bar{X} - \mu_0)).
$$

### 最尤推定量を用いる尤度比検定

上と同じ正規分布のモデルを考える．ただし，帰無仮説 $H_0 : \mu = \mu_0$，対立仮説 $H_1 : \mu \neq \mu_0$ の統計的仮説検定を考える．

$\hat \mu$ を対立仮説の下での最尤推定量 $\hat \mu = \bar{X}$ とする．尤度比の対数の2倍 $\lambda = 2 \log L(\hat{\mu}_1, \mu_0)$ は帰無仮説の下で自由度 1 のカイ二乗分布に従う．

P 値は帰無仮説を棄却できる最小の有意水準であるので，自由度 1 のカイ二乗分布の分布関数 $\chi^2(x)$ を用いて，次で得られる．

$$
\text{P-value} = 1-\chi^2(\lambda).
$$

## 数値計算してみる

以下では上の「ネイマン・ピアソンの補題により最強力検定になる尤度比検定」を method1，「最尤推定量を用いる尤度比検定」を method2 と略記する．

R を用いて次のようにシミュレーション用関数を書いた．

```r
pvsimfun <- function(n, mu, mu0){
  pv1 <- numeric(10000)
  pv2 <- numeric(10000)
  for(i in 1:10000){
    x <- rnorm(n, mu)
    xbar <- mean(x)
    pv1[i] <- pnorm(sqrt(n)*(xbar-mu0), lower.tail=FALSE)
    lr2 <- 2*(sum(dnorm(x, xbar, 1, log=TRUE)) - sum(dnorm(x, mu0, 1, log=TRUE)))
    pv2[i] <- pchisq(lr2, df=1, lower.tail = FALSE)
  }
  data.frame(method1=pv1, method2=pv2)
}

```

まずアルファエラーのシミュレーション結果をP値の分布で示す．

![](/images/two_diff_lrtest/pv_alpha.png)

45度の対角線上に乗っていることから，どちらの手法でも P 値が一様分布になっており，名目上の有意水準が保たれていることがわかる．

次に検出力のシミュレーション結果をP値の分布で示す．帰無仮説は $\mu_0 = 0$ で固定した．

![](/images/two_diff_lrtest/pv_beta.png)

method1 のほうが検出力が高いことがわかる．

作図も含めた R のコード全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/two_diff_lrtest.R