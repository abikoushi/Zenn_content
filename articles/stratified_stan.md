---
title: "Stan ユーザーのための層別解析入門の準備"
emoji: "😸"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, Stan]
published: false
---


## まえおき

因果推論の手法として知られる層別解析を勉強する過程で，私自身がつまずいた前提みたいな部分について書いていく．私自身がつまずいた前提みたいな部分について書いていくので，タイトルは「層別解析入門」ではなく「層別解析入門の準備」とした．

結論を先取りすると，おおむね下の箇条書きのようなことを述べる．

- 介入後の結果変数は観察データのモデルの結果変数とは（パラメータを共有する）別の確率変数だと思うのがいい
- 観察データのモデルが決まっても，それだけでは介入後を表すモデルは一意に決まらない
- 因果グラフも加えて仮定するとそれが一意に決まる


## 確率的グラフィカルモデルを作る

まず（疑似）乱数を使ったシミュレーションをやってみたい．次の表で表せるような同時分布からの乱数をR言語（または適当な統計ソフト）を使って得るにはどうすればいいだろうか．

| $Y$| $X$| $Z$|  確率|
|--:|--:|--:|-----:|
|  0|  0|  0| 0.200|
|  1|  0|  0| 0.050|
|  0|  1|  0| 0.200|
|  1|  1|  0| 0.200|
|  0|  0|  1| 0.050|
|  1|  0|  1| 0.245|
|  0|  1|  1| 0.005|
|  1|  1|  1| 0.050|

三つ組の確率変数 $Y$, $X$, $Z$ はそれぞれ5日後時点での回復（0が未回復で1が回復），風邪薬の投与，重症度を表すことにする．

このモンテカルロ・シミュレーションでは（単変量の）2項分布に従う乱数を生成する `rbinom` 関数を使うことにしよう．

条件付き確率の定義から，同時分布は次のようにも書ける．

$$
\begin{aligned}
P(X,Y,Z) = P(Y|X,Z)P(X,Z) \\
= P(Y|X,Z)P(X|Z)P(Z).\tag{1}
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

向きがあり，元の有向非巡回グラフ（Directed Acyclic Graph; DAG）と呼ぶ．

![](/images/stratified_stan/case1_1.jpg)

そして，このような形でモデルが書けたら，事後分布に従う乱数をサンプリングうるための Stan のコードも書くことができる．

```stan
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

さて，(1)のような分解は一通りに限らない．たとえば次のようにしてもいい．

$$
\begin{aligned}
P(X,Y,Z) &= P(Y|X,Z)P(X,Z) \\
&= P(Y|X,Z)P(Z|X)P(X)).\tag{2}
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

![](/images/stratified_stan/case2_1.jpg)

```stan
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

である．`generated quantities` ブロックを次のように書けばよい．

```stan
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
p(x|z)  =\frac{p(z|x)p(x)}{p(z)}
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
generated quantities{
  real<lower=0, upper=1> gamma;
  vector<lower=0, upper=1>[2] psi;
  gamma = phi[1]*(1-delta) + phi[2]*delta;
  psi[1] = (1-phi[2])*delta/(1-gamma);
  psi[2] = phi[2]*delta/gamma;
}
```

## 因果的グラフィカルモデルを作る

ここまでで同時分布は推定できた．仮にこのような分析を行い「この変数 $x$ が 1 のとき結果 $y$ も 1 がでやすくなります」と報告したとしよう．そうしたら次のステップで「じゃあ $x=1$ にしてみよう！」となるのは自然な発想だと思う．

一方でこれまでの考察で同時分布が定まっても $x=1$ に固定したときの分布は一意でないことがわかる．

つまり (1) の方針でサンプリングする場合， (2) の方針でサンプリングする場合で $x=1$ に固定したときの $y$ の分布は異なる．

(1) の方針でのサンプリングを $x=1$ に固定した場合とは次のようなことだ．

```r
  z = rbinom(n, 1, gamma)
  y0 = rbinom(n, 1, xi[cbind(z+1,1)])
  y1 = rbinom(n, 1, xi[cbind(z+1,2)])
```

(2) の方針でサンプリングする場合は，次のようになる．

```r
  z0 = rbinom(n, 1, phi[1])
  z1 = rbinom(n, 1, phi[2])
  y0 = rbinom(n, 1, xi[cbind(z0+1,1)])
  y1 = rbinom(n, 1, xi[cbind(z1+1,2)])  
```

「Xに効果があることがわかりました」という報告が「ただしここでいう効果はXを実際にやったときの効果ではありません」という意味だとすると分析結果の使いみちに困る感じになる．

しかし順序がつけられると一意に定まる．

![](/images/stratified_stan/case1_2.jpg)



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

![](/images/stratified_stan/case1_3.jpg)

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

![](/images/stratified_stan/case2_2.jpg)

潜在結果変数と呼ぶ．

最尤法を用いてもいい．さらにいえば繰り返し計算も必要なく，


## 感想など：統計的因果推論の必要性について

ここまでで述べたように因果グラフを与えた上でないと因果推論ができないのは物足りなく感じるかもしれない（そう感じるのは僕だけかもしれない）．しかし，上の因果効果の符号が逆とか，大きさが小さすぎとかの場面を想像すると，質的に因果関係がわかったとしても，量的な因果関係を知ることは重要と思えるのではないか．

また，因果グラフについて意見がわかれるようなときも，仮定を明示したほうが「この部分がちょっと違った場合は……」のような議論がクリアになるだろう．ついでに，「ちょっと違った」場合どの程度分析結果が変わりうるか，影響を量的に調べる感度分析というのも必要性があるだろう．


## 参考にした文献など

