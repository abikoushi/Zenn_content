---
title: "Stan ユーザーのための操作変数法"
emoji: "🎻"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, Stan, 統計, 因果推論]
published: true
---

## はじめに

いま，$X$ が $Y$ に与える効果を知りたい．次の図のような因果ダイアグラムを考える．

![](/images/instrument_stan/graph1.png)
*未観測の変数を丸で囲んだ因果ダイアグラム*

知りたい効果とは $X$ に介入を行ったときの効果，すなわち $X$ に入るすべてのパスを削除し，$X$ を適当な値に固定したときの効果である．因果推論の大物であるパールは，このような操作ができることを「因果推論の第一法則」と呼んでいる（[著：ジューディア・パール，ダナ・マッケンジー，訳：夏目大『因果推論の科学』文藝春秋社](https://books.bunshun.jp/ud/book/num/9784163915968) などを参照）．

![](/images/instrument_stan/graph3.png)
*介入後の因果ダイアグラム*

因果推論では，介入前のモデルから介入後のモデルを推定するために，両方に共通するパラメータを導入し，これを推定する．

![](/images/instrument_stan/combine1.jpeg)
*因果推論でやろうとすることの模式図*

具体的なパラメータ化として，正規分布 $\mathcal{N}(\mu, \sigma)$ を用いて次のようなものを考えよう．

$$
\begin{aligned}
Y|X,U &\sim \mathcal{N}(\alpha_{yx}X + \alpha_{yu}U,1)\\
X|U &\sim \mathcal{N}(\alpha_{xu}U,1)\\
U &\sim \mathcal{N}(0,1)\\
\end{aligned}
$$

共通パラメータも含めて模式図を描き直すと次のようになる．

![](/images/instrument_stan/combine2.jpeg)
*因果推論でやろうとすることの模式図．共通パラメータも含めて描いた．*

しかし，このモデルは識別不能である．実際，観測可能である $X$ と $Y$ の共分散は次のようになり，$\alpha_{yu}\alpha_{xu}$ という未知パラメータどうしの掛け算が出てくる．

$$
\begin{aligned}
\mathrm{Cov}(X,Y) &= E[XY]-E[X]E[Y]\\
&= \alpha_{yx}E[X^2] + \alpha_{yu} E[XU]\\
&=\alpha_{yx}+\alpha_{yu}\alpha_{xu}
\end{aligned}
$$

しかし，次のように $U$ から独立で $X$ のみに影響を与える都合のいい変数 $Z$ が観測されていたとすると識別可能になる．

![](/images/instrument_stan/graph2.png)
*操作変数法で扱う状況の因果ダイアグラム*

## Stan による操作変数法

$X$, $Y$ に加えて $Z$ も観測されているとして，次のようなモデルを考えよう．$Z$ は平均0, 分散1になるよう正規化されているとする．

$$
\begin{aligned}
Y|X,U &\sim \mathcal{N}(\alpha_{yx}X + \alpha_{yu}U,1)\\
X|Z,U &\sim \mathcal{N}(\alpha_{xz}Z + \alpha_{xu}U,1)\\
Z & \sim \mathcal{N}(0,1)\\
U & \sim \mathcal{N}(0,1)\\
\end{aligned}
$$

まずは介入前のモデルで乱数をつくってみよう．

```r
rand_pre = function(n, alpha_xu, alpha_xz, alpha_yx, alpha_yu){
  U <- rnorm(n)
  Z <- rnorm(n)
  X <- rnorm(n, alpha_xz*Z+alpha_xu*U)
  Y <- rnorm(n, alpha_yx*X+alpha_yu*U)
  return(data.frame(Y,Z,X))
}

alpha_xu = 0.5
alpha_xz = 0.9
alpha_yx = 0.5
alpha_yu = 0.2

set.seed(1234)
dat <- rand_pre(500, alpha_xu, alpha_xz, alpha_yx, alpha_yu)
```

そして，モデルに現るすべての変数の同時分布を推定するための Stan のコードをいきなり書いてみる．

```stan
data {
  int<lower=0> N;
  vector[N] Y;
  vector[N] X;
  vector[N] Z;
}
parameters {
  vector[N] U;
  real alpha_xz;
  real alpha_yx;
  real<lower=0> alpha_xu;
  real<lower=0> alpha_yu;
}
model {
  for(i in 1:N){
    Z[i] ~ normal(0,1);
    U[i] ~ normal(0,1);
    X[i] ~ normal(alpha_xz*Z[i] + alpha_xu*U[i], 1);
    Y[i] ~ normal(alpha_yx*X[i] + alpha_yu*U[i], 1);
  }
  alpha_xz ~ normal(0, 10);
  alpha_xu ~ normal(0, 10);
  alpha_yx ~ normal(0, 10);
  alpha_yu ~ normal(0, 10);
}
```

ただし，ともに未観測の `alpha_xu`, `alpha_yu` と `U` の掛け算については符号の分の自由度があるから，`alpha_xu`, `alpha_yu` については正の範囲に制約した．つまり，プラス×プラスでもマイナス×マイナスでもプラスになってしまうから，複数のチェインが同じ値に収束するよう一方に制約した．また，`Z[i] ~ normal(0,1);` の行は潜在変数が含まれないので本来は必要ないが，「とにかく全部の同時分布を推定してみる」というニュアンスを出すために入れてある．

介入後のモデルは $X$ に入るすべてのパスを削除し，$X$ を適当な値に固定したときの効果であった．それをシミュレートするために，`generated quantities` ブロックを次のように書いておく．


```stan
generated quantities {
  real TE;
  {
  real Y1;
  real Y0;
  real Uast = normal_rng(0, 1);
  Y1 = normal_rng(alpha_yx+alpha_yu*Uast, 1); // X=1のとき
  Y0 = normal_rng(alpha_yu*Uast, 1);// X=0のとき
  TE = Y1-Y0; // Xを1単位変化させたときYがどの程度変わるか
  }
}
```

`TE` が知りたい効果の事後予測分布である．乱数の作り方から，`TE` の平均は `alpha_yx` になるはずである．

Stan のコードを実行してみると，次のように事後分布に収束していそうな様子が伺える．

![](/images/instrument_stan/trace.png)
*MCMC系列のトレースプロット*

`TE` の事後予測分布のヒストグラムとシミュレーションで設定した真の $\alpha_{yx}$ の値 0.5 を比較する．

![](/images/instrument_stan/hist1.png)
*介入後の効果の事後予測分布．三角形のマーカーは事後予測分布の平均，グレーの点線が真の平均．*

近い値がもとまっていそうなことがわかる．

## オーソドックスな操作変数法

改めて，手計算できるところはもう少し手計算してみる．

$$
\begin{aligned}
\mathrm{Cov}(Y,Z) &= E[YZ]-E[Y]E[Z]\\
&= \alpha_{yx}E[YZ] - E[Y]E[Z]\\
&=\alpha_{yx}\mathrm{Cov}(X,Z)
\end{aligned}
$$

$\mathrm{Cov}(X,Z)$ が得られれば，次式のように $\alpha_{yx}$ の推定量を作ることができる．

$$
\hat{\alpha}_{yx} = \mathrm{Cov}(Y,Z) \, \mathrm{Cov}(X,Z)^{-1}
$$

$Z$ を操作変数と呼ぶ．これが普通の（？）操作変数法として [宮川雅巳『統計的因果推論』（朝倉書店）](https://www.asakura.co.jp/detail.php?book_code=12781&srsltid=AfmBOorK6dX3-enEcoTBr6w80SWGET117fXjczaP5YwvOgNNS3cqAozH)　などの教科書で解説されている．この推定値も先のヒストグラムに重ねてみよう．

![](/images/instrument_stan/hist2.png)
*介入後の効果の事後予測分布．三角形のマーカーは事後予測分布の平均，オレンジの点線が操作変数法による点推定値．*

やはり近い値がもとまっていることがわかる．操作変数法は今回 Stan で行ったようなことの計算コストを劇的に下げる方法だと言えるだろう．

## おわりに

因果推論の手法を最初に学ぼうとしたときは，まずなにをやっているかがわからなかった．介入前後両方のグラフィカルモデルを書き，両者を対応させるというスタイルの下記の文献を見てはじめてなにをやっているかがわかった気がした．

- [Finnian Lattimore & David Rohde - Causal inference with Bayes rule (arXiv)](https://arxiv.org/abs/1910.01510)
- [Finnian Lattimore & David Rohde - Replacing the do-calculus with Bayes rule (arXiv)](https://arxiv.org/abs/1906.07125)
- 著者自身による上の2つ論文の解説記事：[Finnian Lattimore - Causal Inference with Bayes Rule](https://medium.com/gradient-institute/causal-inference-with-bayes-rule-eed8ae45fb2e)

この記事でもそれを踏襲してみたが，どうだろうか．

ところで，今回でいう $Z$ にあたる操作変数は，実際にはそうそう見つからないのではないかという疑問があるかもしれない．
操作変数の具体例としては，[濱谷陸太『極論で語る予防医療』（丸善出版）](https://www.maruzen-publishing.co.jp/book/b10123079.html) 8章で，次のように遺伝子の変異を利用する方法が紹介されている．

- ALDH2の変異があるとお酒が飲めない（暴露因子と関連する）
- 結果変数（心筋梗塞）に曝露因子（飲酒）のみを介して関連する
- ALDH2の変異はランダムで，ALDH2の変異に影響を与えるような因子がない

飲酒の効果を調べるための操作変数として，ALDH2という遺伝子の変異が使えそうだということである．一方で遺伝子変異どうしに関連があって（linkage disequilibrium; 連鎖不平衡というらしい）ALDH2と同時に変異が起こりやすい遺伝子が心筋梗塞に関連しているとか，先祖の遺伝子との関連はあるから地域が交絡因子になるときなど，仮定が成り立たないシナリオもあわせて紹介されている．

また，操作変数と $X$ の相関が小さすぎるときはうまくいかないことにも注意が必要である．これについてもいずれシミュレーションをやってみたいが，今回はここまでで終わりにする．


使用した R のコード全体は以下に置く：

https://github.com/abikoushi/Zenn_content/blob/main/R/instrument_stan.R

