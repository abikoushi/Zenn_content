---
title: "Stan ユーザーのための層別解析入門の準備"
emoji: "💉"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, Stan]
published: false
---

## はじめに

因果推論の手法として知られる層別解析を勉強する過程で，私自身がつまずいた前提みたいな部分について書いていく．前提みたいな部分というのは，そもそもこれなんの計算してるの？　というようなことで，層別解析のやり方についての解説まで至らなかったので，タイトルは「層別解析入門」ではなく「層別解析入門の準備」とした．

結論を先取りすると，おおむね下の箇条書きのようなことを述べる．

- 介入後の結果変数は観察データのモデルの結果変数とは（パラメータを共有する）別の確率変数だと思うのがいい
- 観察データのモデルが決まっても，それだけでは介入後を表すモデルは一意に決まらない
- 因果グラフも加えて仮定するとそれが一意に決まる


確率的グラフィカルモデルを用いた解説のスタイルは主に下記の文献を参考にした．

- [Finnian Lattimore & David Rohde - Causal inference with Bayes rule (arXiv)](https://arxiv.org/abs/1910.01510)
- [Finnian Lattimore & David Rohde - Replacing the do-calculus with Bayes rule (arXiv)](https://arxiv.org/abs/1906.07125)
- [Finnian Lattimore - Causal Inference with Bayes Rule（著者自身による上の2つ論文の解説記事）](https://medium.com/gradient-institute/causal-inference-with-bayes-rule-eed8ae45fb2e)


## 確率的グラフィカルモデルを作る

まず（疑似）乱数を使ったシミュレーションをやってみたい．次の表で表せるような同時分布からの乱数をR言語（または適当な統計ソフト）を使って得るにはどうすればいいだろうか．

| $y$| $x$| $z$|  確率|
|--:|--:|--:|-----:|
|  0|  0|  0| 0.200|
|  1|  0|  0| 0.050|
|  0|  1|  0| 0.200|
|  1|  1|  0| 0.200|
|  0|  0|  1| 0.050|
|  1|  0|  1| 0.245|
|  0|  1|  1| 0.005|
|  1|  1|  1| 0.050|

三つ組の確率変数 $y$, $x$, $z$ はそれぞれ単位時間経過時点での症状の回復（0が未回復で1が回復），薬の投与（0が薬なしで1が薬あり），重症度（0が軽症で1が重症）を表すことにする．

このモンテカルロ・シミュレーションでは（単変量の）2項分布に従う乱数を生成する `rbinom` 関数を使うことにしよう．

条件付き確率の定義から，同時分布は次のようにも書ける．

$$
\begin{aligned}
p(x,y,z) = p(y|x,z)p(x,z) \\
= p(y|x,z)p(x|z)p(z).\tag{1}
\end{aligned}
$$

右辺は単変量の2項分布の積で表せた．2項分布のパラメータを下記のリストのようにおく．

- $p(y=1|z,x) = \xi_{z,x}$
- $p(x=1|z) = \psi_x$
- $p(z=1) = \gamma$ 

これは $(Y, X, Z)$ が単変量の `rbinom` 関数でサンプリングできることを意味する．R のコードはたとえば次のようになる．

```r
rand_case1 = function(n, xi, psi, gamma){
  ## draw x given by z
  z = rbinom(n, 1, gamma)
  x = rbinom(n, 1, psi[z+1]) #R の配列のインデックスは 1 からはじまるため，1 足して z=0 のときを 1行目, z=1 のときを 2 行目にしている
  y = rbinom(n, 1, xi[cbind(z+1,x+1)])
  return(data.frame(Y=y, X=x, Z=z))
}
```
$z$ の実現値が得られれば $x$ もサンプルできるという従属関係を $z \rightarrow x$ と表すことにすると，同時分布は下の図のようなグラフィカルモデルでも表せる．

![](/images/stratified_stan/case1_1.jpg)

このように向きがあり，元のノードに戻るパスがないグラフを有向非巡回グラフ（Directed Acyclic Graph; DAG）と呼ぶ．

そして，このような形でモデルが書けたら，事後分布に従う乱数をサンプリングうるための Stan のコードも書くことができる．

```stan
#case1.stan
data{
  int N;
  array[N] int<lower=0,upper=1> Y;
  array[N] int<lower=0,upper=1> Z;
  array[N] int<lower=0,upper=1> X;
  real<lower=0> alpha;
}
parameters{
    matrix<lower=0, upper=1>[2,2] Xi;
    vector<lower=0, upper=1>[2] psi;
    real<lower=0, upper=1> gamma;
}
model{
  for(i in 1:N){
    Z[i] ~ bernoulli(gamma);
    X[i] ~ bernoulli(psi[Z[i]+1]);
    Y[i] ~ bernoulli(Xi[Z[i]+1, X[i]+1]);
  }
  to_vector(Xi) ~ beta(alpha, alpha);
  gamma ~ beta(alpha, alpha);
  psi ~ beta(alpha, alpha);
}
```

`alpha` は事前分布のパラメータとして導入したが今回はあまり活躍せず，ずっと 1 にしている．

さて，(1)のような分解は一通りに限らない．たとえば次のようにしてもいい．

$$
\begin{aligned}
p(x,y,z) &= p(y|x,z)P(x,z) \\
&= p(y|x,z)p(z|x)P(x)).\tag{2}
\end{aligned}
$$

ここではパラメータを下記のリストのようにおく．

- $p(y=1|z,x) = \xi_{z,x}$
- $p(z=1|x) = \phi_x$
- $p(x=1) = \delta$ 

これに対応する R のコードはたとえば次のように書ける．

```r
rand_case2 = function(n, xi, phi, delta){
  ## draw z given by x
  x = rbinom(n, 1, delta)
  z = rbinom(n, 1, phi[x+1])
  y = rbinom(n, 1, xi[cbind(z+1,x+1)])  
  return(data.frame(Y=y,X=x,Z=z))
}
```

グラフィカルモデルは下の図のようになる．

![](/images/stratified_stan/case2_1.jpg)

Stan のコードは次のように書ける．

```stan
#case2.stan
data{
  int N;
  array[N] int<lower=0,upper=1> Y;
  array[N] int<lower=0,upper=1> Z;
  array[N] int<lower=0,upper=1> X;
  real<lower=0> alpha;
}
parameters{
    matrix<lower=0, upper=1>[2,2] Xi;
    vector<lower=0, upper=1>[2] phi;
    real<lower=0, upper=1> delta;
}
model{
  for(i in 1:N){
    Z[i] ~ bernoulli(phi[X[i]+1]);
    X[i] ~ bernoulli(delta);
    Y[i] ~ bernoulli(Xi[Z[i]+1, X[i]+1]);
  }
  to_vector(Xi) ~ beta(alpha, alpha);
  delta ~ beta(alpha, alpha);
  phi ~ beta(alpha, alpha);
}
```

(1) と (2) は同じ同時分布を表している．変数変換によって $\psi$ と $\gamma$ から $\phi$ と $\delta$ をつくることもできる．

条件付き確率の定義から $p(x)$ は次のようにも表せる．

$$
\begin{aligned}
p(x) &= \sum_z p(x,z) \\
&= \sum_z p(x|z)p(z)
\end{aligned}
$$

この関係から $\delta$ は次式のように表せる．

$$
\begin{aligned}
\delta &= p(x=1) \\
&= \psi_0 \cdot (1-\gamma)+\psi_1 \cdot\gamma.
\end{aligned}
$$

また， $p(z|x)$ は次のように表せる．

$$
\begin{aligned}
p(z|x) &= \frac{p(x,z)}{p(z)}\\
& =\frac{p(x|z)p(x)}{p(z)}
\end{aligned}
$$

この関係を用いて $\phi$ は

$$
\begin{aligned}
\phi_0 &= p(z=1|x=0) \\
&= \frac{(1-\psi_1)\gamma}{1-\delta}
\end{aligned}
$$

および，

$$
\begin{aligned}
\phi_1 &= p(z=1|x=1) \\
&= \frac{\psi_1 \gamma}{\delta}
\end{aligned}
$$

である． $\phi$ と $\delta$ についても事後分布に従がうサンプルを得るには　`generated quantities` ブロックを次のように書けばよい．

```stan
#case1.stan
generated quantities{
  real<lower=0, upper=1> delta;
  vector<lower=0, upper=1>[2] phi;
  delta = psi[1]*(1-gamma) + psi[2]*gamma;
  phi[1] = (1-psi[2])*gamma/(1-delta);
  phi[2] = psi[2]*gamma/delta;
}
```

同様に $\phi$ と $\delta$ から $\psi$ と $\gamma$ をつくることもできる．

$p(z)$ は次のように表せる．

$$
p(z) = \sum_x p(z|x)p(x)
$$

これより $\gamma$ は，

$$
\begin{aligned}
\gamma &= p(z=1) \\
&= \phi_0 \cdot (1-\delta)+\phi_1 \cdot\delta.
\end{aligned}
$$

また $p(x|z)$ は次のように表せる．

$$
p(x|z) = \frac{p(z|x)p(x)}{p(z)}
$$

この関係を用いて $\psi$ は

$$
\begin{aligned}
\psi_0 &= p(x=1|z=0) \\
&= \frac{(1-\phi_0)\delta}{\gamma}
\end{aligned}
$$

および，

$$
\begin{aligned}
\psi_1 &= p(x=1|z=1) \\
&= \frac{(1-\phi_1)\delta}{\gamma}
\end{aligned}
$$

である．

```r
#case2.stan
generated quantities{
  real<lower=0, upper=1> gamma;
  vector<lower=0, upper=1>[2] psi;
  gamma = phi[1]*(1-delta) + phi[2]*delta;
  psi[1] = (1-phi[2])*delta/(1-gamma);
  psi[2] = phi[2]*delta/gamma;
}
```

次の図からMCMCの系列が定常分布に収束した様子がうかがえる．

![](/images/stratified_stan/traceplot.png)

はほぼ同じ事後分布が得られている．

![](/images/stratified_stan/param_ecdf.png)

重なっていてほぼ違いが見えない．

|variable | model1.mean| model1.sd| model2.mean| model2.sd| true|
|:--------|-----------:|---------:|-----------:|---------:|----:|
|Xi[1,1]  |        0.22|      0.04|        0.22|      0.04| 0.20|
|Xi[2,1]  |        0.86|      0.03|        0.86|      0.03| 0.83|
|Xi[1,2]  |        0.45|      0.04|        0.45|      0.04| 0.50|
|Xi[2,2]  |        0.88|      0.06|        0.88|      0.06| 0.91|
|psi[1]   |        0.58|      0.03|        0.58|      0.03| 0.62|
|psi[2]   |        0.18|      0.03|        0.18|      0.03| 0.16|
|phi[1]   |        0.51|      0.03|        0.51|      0.03| 0.54|
|phi[2]   |        0.14|      0.02|        0.14|      0.02| 0.12|
|delta    |        0.44|      0.02|        0.44|      0.02| 0.35|
|gamma    |        0.34|      0.02|        0.35|      0.02| 0.46|


## 因果的グラフィカルモデルを作る

ここまでで同時分布は推定できた．仮にこのような分析を行い「この変数 $x$ が 1 のとき結果 $y$ も 1 がでやすくなります」と報告したとしよう．そうしたら次のステップで「じゃあ $x=1$ にしてみよう！」となるのは自然な発想だと思う．

一方でこれまでの考察で同時分布が定まっても $x=1$ に固定したときの分布は一意でないことがわかる．

つまり (1) の方針でサンプリングする場合， (2) の方針でサンプリングする場合で $x=1$ に固定したときの $y$ の分布は異なる．

(1) の方針でのサンプリングを $x=1$ に固定した場合とは次のようなことだ．

```stan
#case1
generated quantities{
  int D;
  {
    int Zast = bernoulli_rng(gamma)+1;
    int Yast0 = bernoulli_rng(Xi[Zast,1]);
    int Yast1 = bernoulli_rng(Xi[Zast,2]);
    D = Yast1 - Yast0;
  }
}
```

![](/images/stratified_stan/case1_2.jpg)

(2) の方針でサンプリングする場合は，次のようになる．

```stan
#case2
generated quantities{
int D;
{
  int Yast0;
  int Yast1;
  int Zast0;
  int Zast1;
  Zast0 = bernoulli_rng(phi[1])+1;
  Zast1 = bernoulli_rng(phi[2])+1;
  Yast0 = bernoulli_rng(Xi[Zast0,1]);
  Yast1 = bernoulli_rng(Xi[Zast1,2]);
  D = Yast1 - Yast0;
}
}
```

下の図のような

![](/images/stratified_stan/case1_3.jpg)

「Xに効果があることがわかりました」という報告が「ただしここでいう効果はXを実際にやったときの効果ではありません」という意味だとすると，分析結果の使いみちに困る感じになる．

しかし順序がつけられると一意に定まる．仮に


仮に次のように因果グラフが与えられたとき，

![](/images/stratified_stan/case2_3.jpg)

アルゴリズミックに得られる．

![](/images/stratified_stan/case2_2.jpg)

$y^{\ast}$ を潜在結果変数（potential outcome）と呼ぶ．

|variable    |  mean|   sd|
|:-----------|-----:|----:|
|D1 (model1) |  0.16| 0.60|
|D2 (model1) | -0.03| 0.71|
|D1 (model2) |  0.15| 0.61|
|D2 (model2) | -0.04| 0.71|

 0.223

 0.008

実はこの設定ではMCMCを用いる必要はなく，事後分布はベータ分布であり，事後分布のパラメータも閉じた形で求まる．もちろん最尤法を用いてもいい．

最尤法を用いた場合の具体的な計算の手順は例えば，

- [佐藤俊哉・松山裕，交絡という不思議な現象と交絡を取りのぞく解析—標準化と周辺構造モデル—，計量生物学 / 32 巻 (2011) ](https://www.jstage.jst.go.jp/article/jjb/32/Special_Issue/32_Special_Issue_S35/_article/-char/ja/) 
- [佐藤俊哉『宇宙怪人しまりす統計よりも重要なことを学ぶ』（朝倉書店）](https://www.asakura.co.jp/detail.php?book_code=12297&srsltid=AfmBOooBWqDbRmCFU-gQMF1sM3If4E48Pir3moPOL8H3hICqFBQogt6v)の第4話
- [林岳彦『はじめての統計的因果推論』（岩波書店）](https://www.iwanami.co.jp/book/b639904.html) の4.1節「層別化と標準化で揃える」

などに解説がある．

ちなみに，『宇宙怪人しまりす統計よりも重要なことを学ぶ』では「層別解析」だとサブグループ解析などと混同しやすいことから，「層化統合解析」という用語を提案している．たしかにそのとおりだと思うので今後は「層化統合解析」の方を採用するかもしれない．今回はこの文章だけで完結するような書き方ではない （あまり self-contained でない）ので他の文献を参照しやすいよう，よく使われる「層別解析」を採用した．


## 感想など

統計的因果推論の必要性について：

- 最初，因果グラフを与えた上でないと因果推論ができないのは物足りなく感じたことがある．しかし，上の因果効果の符号が逆とか，大きさが小さすぎとかの場面を想像すると，質的に因果関係がわかったとしても，量的な因果関係を知ることは重要と思えた．
- 因果グラフについて意見がわかれるようなときも，仮定を明示したほうが「この部分がちょっと違った場合は……」のような議論がクリアになるだろう．
- ついでに，「ちょっと違った」場合どの程度分析結果が変わりうるか，影響を量的に調べる感度分析というのも必要性があるだろう．


もう少し細かい話：

- 数理統計の教科書では推定量のバイアス（bias）はパラメータ $\theta$ とその推定量 $\hat \theta$ の差 $\theta - E[\hat{\theta}]$ と定義される．因果推論の文脈でバイアスを減らすことが強調される場合は，もう少し広い意味でのバイアス（標本サイズを大きくしても小さくならないようなバイアス）と捉えるべきだろう．つまり不偏推定量を使わなくとも因果推論の目的に照らして矛盾はないだろう．
- また，上で標準誤差を見たような場面でP値を見ることも因果推論の目的とはまったく対立しないだろう．
