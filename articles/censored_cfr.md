---
title: "打ち切りのあるときの致命割合"
emoji: "✂️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学]
published: false
---

## はじめに

このノートは [西浦博　編著『感染症流行を読み解く数理』（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) の第8章（担当：西浦博）を参考にしたものです．

致命割合（case fatality ratio; CFR） は 「感染者のうち， 感染というイベントを通じて死亡した者の割合」 を指す．

これは $t$ 時点までの累積の感染者数 $C_t$ と累積の感染者数中の死亡者数 $D_t$ の比 $D_t/C_t$ で簡単に推定できそうに思える．しかし，$t$ 時点までの死亡者数は $t$ 時点までの打ち切りデータであるため，推定量 $D_t/C_t$ は致命割合を過小評価するおそれがある．

以下ではこの問題への対処を考えるため，簡単なシミュレーションを行う．


## シミュレーション

感染というイベントの発生時刻 $t^c_i$ は区間 $[0,1), [1,2), \ldots , [t-1, t)$ のそれぞれで区分的に定数の強度 $\lambda(t)$ を持つ非定常ポアソン過程に従うとする．また，感染から $x_i$ 遅れて死亡というイベントが発生するとする．

すなわち次のようなデータ生成過程を設定し，シミュレーションを行う．

$$
\begin{aligned}
t^c_i &\sim \mathrm{NHPP}(\lambda(t))\\
t^d_i &= t^c_i+x_i\\
(x_i \mid d_i=1) &\sim F(\cdot)\\
d_i &\sim \mathrm{Bernoulli}(p)
\end{aligned}
$$

$t^d_i <t$ を満たす $t^c_i$ の数を $C_t$ ，$t^d_i <t$ を満たす $t^d_i$ の数を $D_t$ とする．

イメージしやすいように，一回のシミュレーション結果を図示してみよう．

```r
Tmax <- 100
set.seed(1234)
n <- rpois(1,10*Tmax) #レート10の定常ポアソン過程
ti_c <- sort(runif(n, 0, Tmax))
d <- rbinom(n, 1, 0.2) #真の致命割合は0.2
ti_d <- sort(ti_c + ifelse(d==1, rlnorm(10000, 2, 0.5), Inf)) #Fは対数正規分布とした

Cases <- sapply(1:Tmax, function(t)sum(ti_c<=t))
Deaths <- sapply(1:Tmax, function(t)sum(ti_d<=t))
```

![](/images/censored_cfr/obs.png)
*Casesが $C_t$，Deaths が$D_t$．Cases に対して Deaths は緩やかな傾きで，遅れて立ち上がってくる．*

さて，さらに $f_t = F(t)-F(t-1)$ とおく． また $\lambda_t$ を次のようにおく．

$$
\lambda_t = \int^{t}_{t-1}\lambda(u) \, du.
$$

このとき $D_t$ は次のパラメータを持つポアソン分布に従う．

$$
\begin{aligned}
D_t \sim \mathrm{Poisson}\left(\sum_{h=0}^{\infty}pf_h\lambda_{t-h}\right).
\end{aligned}
$$

$t-1 \le t^c_i <t$ を満たす $t^c_i$ の数を $c_t$ とする．$\lambda_t$ の推定量を $c_t$ そのものとして $A_t = \sum_{h=0}^{\infty}f_h c_{t-h}$ と置くと，$p$ についての対数尤度は次のように書ける．

$$
L(p) = D_t \log(pA_t) - (pA_t) + \mathrm{Const.}
$$

最尤推定量は次のように得られる．

$$
\hat{p} = D_t/A_t.
$$

これは [西浦博　編著『感染症流行を読み解く数理』（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) の 8.4 式と一致する．ただし [『感染症流行を読み解く数理』](https://www.nippyo.co.jp/shop/book/8827.html) では上述のようなモデルは明記されていない．ここではシミュレーションのため，私の判断で補った．

また[『感染症流行を読み解く数理』](https://www.nippyo.co.jp/shop/book/8827.html) では $\hat{p}$ に特別な名前はついていないが，ここでは区別のため aCFR と呼ぶことにする．一方，$D_t/C_t$ を cCFR と呼ぶことにする．

aCFR と cCFR をシミュレーションで比較しよう．

![](/images/censored_cfr/fit.png)
*x軸の時点までで評価した致命割合の比較．100回のシミュレーション結果．*

aCFR は真の致命割合を中心に分布していることがわかる．ただし常に保守的な推定値を与え，最初期は1を超える値が推定されることもある．

R のコード全体は以下：

https://github.com/abikoushi/Zenn_content/blob/main/R/censored_cfr.R