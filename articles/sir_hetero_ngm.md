---
title: "人口の異質性を加味した SIR モデルに基づく次世代行列の構築 "
emoji: "👨‍👦‍👦"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 微分方程式, 線形代数]
published: false
---

## このノートについて

[西浦博　編著『感染症流行を読み解く数理』（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) の第１章の「人口の異質性を加味した SIR モデル」の節を，自分なりに数値例を補いながら読んでいくものです．この第１章の筆者は小林鉄郎，西浦博です．ただし，（当たり前かもしれませんが，）この文書の責任は私にあります．本文中のコードは R 言語です．

## 人口の異質性を加味した SIR モデル

人口に $n$ 個のサブグループがあるとして，$a=1,\ldots, n$ について次の連立微分方程式で表されるモデルを考える．

$$
\begin{aligned}
\frac{d}{dt}S_a(t) &= -S_a(t) \sum_{b=1}^n \beta_{ab} I_b(t)\\
\frac{d}{dt}I_a(t) &= S_a(t) \left(\sum_{b=1}^n \beta_{ab} I_b(t) \right) - \gamma I_b(t)\\
\frac{d}{dt}R_a(t) &= \gamma I_a(t)
\end{aligned} \tag{1}
$$

これは人口の異質性を加味した SIR モデルであり，SIRモデルには感受性者 $S$，感染期 $I$，回復・免疫保持 $R$ の3つの区画がある．各パラメータは次のような意味を持つ．

- $\beta_{ab}$: サブグループ $b$ の感染期の個体によるサブグループ $a$ の個体への感染伝播率
- $\gamma$: 回復率（$1/\gamma$ が平均感染期間）

今回は $S_a(t)+I_a(t)+R_a(t)=1$ になるよう，サブグループごとに正規化されているとする．

 $S(t) = (S_1(t),S_2(t), \ldots, S_n(t))'$, $I(t) = (I_1(t),I_2(t), \ldots, I_n(t))'$, $T(t) = (I_1(t),I_2(t), \ldots, R_n(t))'$， さらに $\beta = (\beta_{ab})$ とまとめて置く．

$$
\begin{aligned}
\frac{d}{dt}S(t) &= -S(t) \circ \beta I(t)\\
\frac{d}{dt}I(t) &= S(t) \circ \beta I(t) - \gamma I(t)\\
\frac{d}{dt}R(t) &= \gamma I(t)
\end{aligned}
$$

ここで $\circ$ はアダマール積（ベクトルの要素ごとの積）とした．

![](/images/sir_hetero_ngm/SIR1.png)

## 次世代行列の導入

$S_a(t) \approx 1$ ($a=1,\ldots g$) が成り立つ感染流行初期を考える．このとき，感染 $I$ についての方程式は次のようになる．

$$
\frac{d}{dt}I(t) = \beta I(t) - \gamma I(t)
$$

これを次のような行列の和の形に分解しよう．

$$
\frac{d}{dt}I(t) = (T + \Sigma)I(t).
$$

ここで $T$ は新規感染の発生に関する項に注目した係数である．

$$
T = \beta.
$$

また $\Sigma$ は新規感染以外の状態の変化（ここでは回復）に関する項のみに注目した係数である．

$$
\Sigma = \begin{pmatrix}
 -\gamma &  0 & \cdots &  0\\
0 &  -\gamma & \cdots &  0\\
\vdots& \vdots  & \ddots &  \vdots\\
0 &  0 & \cdots &  -\gamma\\
\end{pmatrix}
.
$$

初期の感染は時間 $t$ について指数関数 $\exp((T+\Sigma)t)$ で近似できる．この指数関数と先程の SIR モデルの解を重ねてみる．

![](/images/sir_hetero_ngm/SIR_exp.png)

この分解に基づき行列 $K_L$ を次のように定義する．

$$
K_L = -T \Sigma^{-1}.
$$

$K_L$ は「感染者が各状態で過ごす平均時間」と「その状態にいるときに新規感染を発生させる速度」を掛け合わせることで，「1人の感染者が生み出す新規感染者の総数」を表している．基本再生算数は $K_L$ の最大固有値で定義できる．

$K_L$ の要素を書き下すと，

$$
K_L = 
\begin{pmatrix}
\beta_{11}/\gamma &  \beta_{12}/\gamma & \cdots &  \beta_{1n}/\gamma\\
\beta_{21}/\gamma &  \beta_{22}/\gamma & \cdots &  \beta_{2n}/\gamma\\
\vdots&  \vdots & \ddots &  \vdots\\
\beta_{n1}/\gamma &  \beta_{n2}/\gamma & \cdots &  \beta_{nn}/\gamma\\
\end{pmatrix}.
$$

以降では $K_L$ の $(a,b)$ 成分を $R_{ab}$ と置く．

## 最終規模方程式

(1) の $R(t)$ についての方程式より，

$$
I_a(t) = \frac{1}{\gamma }\frac {d}{dt}R_a(t).
$$

これを (1) の $S_a(t)$ についての方程式に代入すると，次の方程式を得る．

$$
\frac{1}{S_a(t)} {\frac {d}{dt}}S_a(t)=-  \sum_{b=1}^n\frac{\beta_{ab}}{\gamma}\frac {d}{dt}R_a(t).
$$

$z_a = \lim_{t \to \infty}R_a(t)$ と置き，$\lim_{t \to \infty}S_a(t) = 1-z_a$ に注意して，両辺を区間 $[0, \infty)$ で積分することで次の方程式を得る．

$$
\log(1-z_a)= \sum_{b=1}^n\frac{\beta_{ab}}{\gamma} z_b.
$$

整理すると，

$$
z_a= 1- \exp\left(\sum_{b=1}^nR_{ab} z_b \right).
$$

この方程式を満たす $z_a$ が最終規模である．

![](/images/sir_hetero_ngm/SIR_R.png)