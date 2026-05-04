---
title: "SIRモデルと基本再生産数"
emoji: "💉"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 微分方程式, 疫学]
published: false
---

## このノートについて

[西浦博編著『感染症流行を読み解く数理』（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) の第１章のはじめの方を自分なりに適当と思われるRによる数値例を補いながら読んでいくものです．主に参照した第１章の筆者は小林鉄郎，西浦博です．ただし，当たり前かもしれませんが，この文書の責任は私にあります．

## SIRモデル

SIRモデルは感染症の流行過程を記述するモデルで，次の連立微分方程式で表される．

$$
 \begin{aligned}
 {\frac {d}{dt}}S(t)&=-\beta S(t)I(t) \\
 {\frac {d}{dt}}I(t)&=\beta S(t)I(t)-\gamma I(t)\\
 {\frac {d}{dt}}R(t)&=\gamma I(t)
 \end{aligned}
 \tag{1}
$$

$S(t)$ は $t$ 時点でこれから感染する可能性のある人の数，$I(t)$ は感染者数，$R(t)$ は除外数——すなわち病気から回復して免疫を獲得した人などもう感染する可能性のない人の数に対応する． $\beta$, $\gamma$ はパラメータで，$\beta$ は感染者 $I(t)$ 一人あたりかつ単位時間あたりの感染率，$\gamma$ は除外率である．今回は $S(t)+I(t)+R(t)=1$ になるよう正規化されているとする．

SIRモデルの解をプロットしてみよう．微分方程式を数値的に解くRのパッケージ `deSolve` がある．SIRモデルを数値的に解くには次のようにすればよい．

```
library(deSolve)

SIRmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS <- - beta*S*I
    dI <- beta*S*I - gamma*I
    dR <- gamma*I 
    return(list(c(dS,  dI, dR)))
  })
}

pars  <- c(beta=1.5, gamma=1)

ini = c(S=0.99, I=0.01, R=0)
ode_out1 <- ode(y=ini, times=times, func=SIRmod, parms=pars)
```

![](/Zenn_content/images/sir_model_and_r0/SIR_sol.png)
*SIRモデルの解 $S(t), I(t), R(t)$ のグラフ．横の点線については後述する．*

パラメータは `pars  <- c(beta=1.5, gamma=1)` の部分で　$\beta=1.5, \gamma=1$ と与えた．

さて，感染症の流行時に感染が広がるのか収束していくのかには特に興味がある．そこで初期値をいくつか振り，感染者数 $I(t)$ をプロットしてみる．

![](/Zenn_content/images/sir_model_and_r0/SIR_sol.png)
*SIRモデルの解 $I(t)$ についてのグラフ．*

ある値は
この「ある値」は $\gamma / \beta$ である．

初期状態での変化量が正になるのは, $I(t)$ についての方程式

$$
 {\frac {d}{dt}}I(t)= (\beta S(t)-\gamma) I(t)
$$

に着目すると，次のときである．

$$
\beta S(0)-\gamma > 0.
$$

すなわち

$$
 S(0) > \gamma/\beta
$$

のときである．これが $\gamma/\beta$ を境に解の様子が大きく変わる理由である．初期状態では $I(0)$ が 0 に近いと考えられるから， $R(0) \approx 1-S(0)$ とすると，次のようにも表せる．

$$
 R(0) < 1 - \gamma/\beta.
$$

SIRモデル内においては，ワクチンで免疫を獲得した人が $ 1 - \gamma/\beta$ より多ければ感染が拡大しないだろうと見積もれることを意味する．


## 基本再生産数

先の不等式はまた，次のように表すこともできる. 

$$
\frac{\beta}{\gamma} S(0) > 1. 
$$

いま，$S(0)$ がほぼ 1，$R(0)$ が0の状況を考え，左辺の $\beta/\gamma$ を基本再生産数と呼ぶ．以下では基本再生産数を $R_0 = \beta/\gamma$ とおく．この $R_0$ は $R(0)$ とは記号が似ているだけで別のものである．

今度は，感染が広がるときは最終的にどのくらいまで広がるのかに興味がある．

(1) の $R(t)$ についての方程式より

$$
I(t) = \frac{1}{\gamma }\frac {d}{dt}R(t)
$$

を (1) の $S(t)$ についての方程式に代入すると，次の方程式を得る．

$$
\frac{1}{S(t)} {\frac {d}{dt}}S(t)=-R_0  \frac {d}{dt}R(t).
$$

$z = \lim_{t \to \infty}R(t)$ と置き，$\lim_{t \to \infty}S(t) = 1-z$ に注意して，両辺を区間 $[0, \infty)$ で積分することで次を得る．

$$
\log (1-z) -\log S(0) = -R_0  (z - R(0)).
$$

整理すると，

$$
R_0 = - \frac{\log (1-z) -\log S(0)}{z-R(0)}. 
$$

実はこれを $z$ について解いた値が最初の図に引いていた点線の意味である．いま $R(0)$ は0，$S(0)$ はほぼ1の状況を考えていたので，上の等式を次のように表すことができる．

$$
R_0 = - \frac{\log (1-z)}{z}.
$$

これは最終規模方程式と呼ばれる．SIRモデル内においては，感染が拡大したときの最終的な規模が，基本再生産数 $R_0$ のみによって決まることを意味する．$R_0$ を動かして $z$ をプロットしてやると次のようになる．

![](/Zenn_content/images/sir_model_and_r0/size_by_R0.png)
*基本再生産数ごとの最終的な感染者数*

## 落穂拾い

微分方程式で表されるモデルの傾きが 0 になる点の集合をヌルクラインと呼ぶ．ヌルクラインを境に解の振る舞いが大きく変わるというのは，SIRモデルに限ったことではない．このことは [Hirsch・Smale・Devaney 力学系入門（共立出版）](https://www.kyoritsu-pub.co.jp/book/b10003811.html) などに詳しい．