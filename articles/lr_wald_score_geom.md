---
title: "尤度比検定, ワルド検定, スコア検定：幾何分布の例"
emoji: "🚌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 仮説検定]
published: false
---

## 概要

[尤度比検定, ワルド検定, スコア検定から定まる信頼区間：二項分布の例](https://zenn.dev/abe2/articles/lr_wald_score_binom)

の続き的な内容．

二項分布の例は教科書にも比較的よく出ているので，あまり出てないであろう幾何分布でもやってみようというもの．

こんな感じで自分でも好きな分布（モデル）で検定や信頼区間を作ってみようかなと思ってもらえたら大変うれしい（……が，そのためには本当は解析的に解けないような場合の例も必要だと思う）．


## 動機付けのための例

あるソーシャルゲーム（ソシャゲ）ではユーザーが1日あたりどのくらいの確率でログインしているか知りたいと思っている．

マーケティングの分野では $t_n$ をリセンシー，$n$ をフリクエンシーとして，リセンシー・フリクエンシー・マネタリーバリュー（購買金額）をあわせて色々考えることをRFM分析と呼ぶ場合がある．

幾何分布という用語にあまり馴染みのない方のため，これまで述べたことを再度図にまとめておく．

ここではなるべくイメージが持ちやすいようにこのような例にしてみたが，以降特別ソシャゲの話題が出てくるわけではないので，医学データに興味のある人はログインの代わりに症状の再発とか，品質管理に興味のある人は部品の故障とか，自分の興味のある題材に読み替えてもらえると嬉しい．


## 地道に手計算パート

尤度に基づく方法を考える.

上の設定のモデルで，尤度関数は

$$
L(p) = \left( \prod_{i=1}^n (1-p)^{t_i-t_{i-1}-1} p \right) ((1-p)^{w - t_n})
$$

と書ける. 推定したい未知パラメータは $p \in [0,1]$ である.

対数尤度関数を少し整理する．

$$
\begin{aligned}
\log L(p) &= \left( \sum_{i=1}^n (t_i-t_{i-1}-1)\log(1-p) + \log p \right) + (w - t_n)\log(1-p) \\
& = \{(t_1-t_{0}-1)+(t_2-t_{1}-1)+ \cdots + (t_n-t_{n-1}-1) + (w - t_n)\}\log(1-p) + n \log p  \\
&= (w-n) \log(1-p) + n \log p. 
\end{aligned}
$$

ここまで求めると，あとの尤度比検定・スコア検定・ワルド検定の導出は2項分布のときとまったく同じ計算になる．

[尤度比検定, ワルド検定, スコア検定から定まる信頼区間：二項分布の例](https://zenn.dev/abe2/articles/lr_wald_score_binom) の結果で $x$ を $n$ に，$n$ を $w$ に，$\theta$ を $p$ に置き換えればよい．

置き換えるだけだが2つの記事を見比べるのが面倒かもしれなので，この記事の最後にまとめて書いておく．

最尤推定量は

$$
\hat p = n/w
$$

と求まる.

改めて考えると，この問題は「$w$ 回の試行のうちログインした回数 $n$ の分布」を扱っているのと同じなので，これはあたりまえのことであった．

また，尤度を計算するのには履歴 $t_i$ をすべてコンピュータ上のメモリに乗せる必要はなく，$w$ と $n$ がわかれば十分である．
このような統計量を一般化すると十分統計量という考え方になる．


## Rによる実装例とシミュレーション

尤度比検定に基づく p 値を返す関数は例えば次のように書ける.

```R
LRtest <- function(p,x,n){
  phat <- x/n
  lp <- 2*(dbinom(x, n, phat ,log = TRUE)-
             dbinom(x, n, p ,log = TRUE))
  pchisq(lp, 1, lower.tail = FALSE)
}
```

Wald 検定に基づく p 値を返す関数は例えば次のように書ける.

```
waldtest <- function(p,x,n){
  phat <- x/n
  se <- sqrt(phat*(1-phat)/n)
  2*pnorm(abs(phat-p)/se, lower.tail = FALSE)
}

pv_w <- sapply(z, waldtest, x=x, n=n)
df_w <- data.frame(p=z, pv=pv_w, method="Wald")

```

`abs` で絶対値を取って片側だけ計算しているので2倍している.

スコア検定に基づく. このことは次のように確かめられる.


```r
set.seed(1234)
res_p <- prop.test(x, n, p = 0.5, correct = FALSE, conf.level = 0.95)

CI_score <- function(x, n, level){
  z <- qnorm(0.5*(1-level), lower.tail = FALSE)
  phat <- x/n
  t_1 <- phat+(z^2)/(2*n)
  t_2 <- z*sqrt(z^2/(4*n^2)+phat*(1-phat)/n)
  c((t_1 - t_2)/(1+(z^2)/n),
    (t_1 + t_2)/(1+(z^2)/n))  
}


print(CI_score(x,n,0.95))
print(res_p$conf.int)
```

```r
> print(CI_score(x,n,0.95))
[1] 0.2581979 0.6579147
> print(res_p$conf.int)
[1] 0.2581979 0.6579147
attr(,"conf.level")
[1] 0.95
```

 p 値を 3 つの検定で比較してみよう.





Rのコード全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/28b54a4fa9641c360d8df0fb285b29e8d5b9f78c/R/lr_wald_score_geom.R




## 尤度比検定, ワルド検定, スコア検定の導出

### 準備

最尤推定量のサンプルサイズ $n$ の大きいときの分布（漸近分布）を知るために，フィッシャー情報量を考えておく.

フィッシャー情報量 $I(\theta)$ は対数尤度関数の2階微分の期待値の符号反転である.

未知パラメータで1階微分した導関数をスコア関数と呼ぶ.

スコア関数 $S(\theta)$ は

$$
S(p) = \frac{d}{dp}\log L(p) = .
$$

フィッシャー情報量は

$$
\begin{align*}
I(\theta) &= -E\left[\frac{d^2}{d\theta^2}\ell(\theta)\right] \\
&= -(E[-X /\theta^2 - (1-X)/(1-\theta)^2])\\
&= -n(-1/\theta - 1/(1-\theta))\\
&= \frac{n}{\theta(1-\theta)}.
\end{align*}
$$

### 尤度比検定

帰無仮説 $\theta=\theta_0$, 対立仮説 $\theta \neq \theta_0$ の尤度比検定では $\theta = \theta_0$ と固定した自由パラメータが0個のモデルと, 自由パラメータが1個のモデルで尤度の比の対数の2倍が自由度 $1-0 = 1$ のカイ二乗分布に従うことを利用する.


検定統計量は

$$
T(X) = 2(\ell(\hat \theta) -\ell(\theta_0))
$$

で, これが自由度1のカイ二乗分布の上側確率 $\alpha$ を与える点以下の範囲の $\theta_0$ が信頼区間である.

しかし閉じた形では求まらなかった.


### Wald 検定


帰無仮説 $\theta=\theta_0$, 対立仮説 $\theta \neq \theta_0$ の Wald 検定は検定統計量

$$
T(X) = |\hat \theta - \theta|\sqrt{I(\hat \theta)}
$$

が帰無仮説の下で漸近的に標準正規分布にしたがうことを利用する. すなわち, 帰無仮説の下で

$$
P_{\theta_0}(T(X) < z_{\alpha/2}) = 1-\alpha
$$

ここで $z_{\alpha/2}$ は標準正規分布の上側確率 $\alpha/2$ を与える分位点である.

信頼区間はこれを $\theta_0$ について解くことで,

$$
C(X) = [C_{-}, C_{+}]
$$

と求まる.

ただし,

$$
C_{-} = \hat \theta - z_{\alpha/2}\sqrt{\hat \theta(1 - \hat\theta)/n},
$$

$$
C_{+} = \hat \theta +z_{\alpha/2}\sqrt{\hat \theta(1- \hat\theta)/n}
$$

とした.


### スコア検定


帰無仮説 $\theta=\theta_0$, 対立仮説 $\theta \neq p_0$ のスコア検定は検定統計量

$$
T(X) = |S(\theta_0)|/\sqrt{I(\theta_0)}
$$

が帰無仮説の下で漸近的に標準正規分布にしたがうことを利用する.

この場合について書き下すと,

$$
T(X) = \sqrt{n} \, \frac{|\hat \theta - \theta_0|}{\sqrt{\theta_0(1- \theta_0)}}.
$$

$T(X) \le z_{\alpha/2}$ となる（受容域に入る）のは,

$$
(\hat \theta - \theta_0)^2 < z_{\alpha/2}^2 \,\frac{\theta_0(1- \theta_0)}{n}
$$

のときである.

この不等式は次のようにも書ける.

$$
\left(1+\frac{z_{\alpha/2}^2}{n}\right)\theta_0^2 - 2\left(\hat\theta\theta_0 + \frac{z_{\alpha/2}}{n}\right)\theta_0 + (\hat \theta)^2  \le  0
$$

なので $\theta_0$ について解くことは2次方程式の解を求めることであり, いくらかの計算の後, 

$$
C(X) = [C_{-} , C_{+}]
$$

ただし,

$$
C_{-} = (1+\frac{z_{\alpha/2}}{n})^{-1} \left( \hat \theta - \frac{z_{\alpha/2}^2}{2n} - z_{\alpha/2}\sqrt{\frac{z_{\alpha/2}}{4n^2} + \frac{\hat \theta (1-\hat \theta)}{n}} \right),
$$

$$
C_{+} = (1+\frac{z_{\alpha/2}}{n})^{-1} \left( \hat \theta + \frac{z_{\alpha/2}^2}{2n} + z_{\alpha/2}\sqrt{\frac{z_{\alpha/2}}{4n^2} + \frac{\hat \theta (1-\hat \theta)}{n}}\right)
$$

とした.
