---
title: "一般分枝過程のマルサス係数（ガンマ分布の場合）"
emoji: "🌿"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, Rcpp, 確率過程]
published: true
---

## はじめに

Nerman, O. (1981). On the convergence of supercritical general (CMJ) branching processes. Zeitschrift für Wahrscheinlichkeitstheorie und verwandte Gebiete, 57(3), 365-395. https://link.springer.com/content/pdf/10.1007/BF00534830.pdf

などを見ると，一般分枝過程（general branching processes）のマルサス係数（Malthusian parameter） $\alpha$ は次の方程式を満たす $\alpha$ として定義されている．

$$
\int_0^\infty e^{-\alpha t} \mu(dt)=1
$$

ここで $\mu(t)$ は時刻 $t$ までに生起するノードの数 $N(t)$ の平均，

$$
\mu(t) = E[N(t)].
$$

これが標準的な定義のようだが，正直これだけだとよくわからないので今回は

1. 最初は始祖となる1ノードからスタート
2. 各ノードは出生後、ガンマ分布に従う待ち時間を繰り返して1人ずつ子を産む（再生過程）

という設定で，具体的にマルサス係数を考えてみる．ここでは特に個体に死亡を仮定せず、無限にノードを産み続けるモデルを考える．

## ガンマ分布の場合のマルサス係数

再生間隔を $X_i$ とし，独立同分布で形状パラメータ $k$ レートパラメータ $\beta$ のガンマ分布に従うとする．

$$
X_i \sim \operatorname{Gamma}(k,\beta)
$$

個体 $j$ についての出生時刻を $S_j$ とする．

$$
S_j=X_1+\cdots+X_j.
$$

今回のモデルでは，

$$
\int_0^\infty e^{-\alpha t} \mu(dt)=
\int_0^\infty e^{-\alpha t} \sum_{n}^{\infty}P(S_n \in dt)
$$

より，マルサス係数は次の方程式を満たす $\alpha$ である．

$$
\sum_{n=1}^{\infty}E[e^{-\alpha S_n}]=1
$$

いま，ガンマ分布のモーメント母関数より，

$$
E[e^{-\alpha X}]=\left(\frac{\beta}{\beta+\alpha}\right)^k.
$$

$q(\alpha)= E[e^{-\alpha X_i}]$ と置くと，独立性から

$$
E[e^{-\alpha S_n}]
= E[e^{-\alpha X_1}]\cdots E[e^{-\alpha X_n}]
=q(\alpha)^n.
$$

したがってマルサス係数を決める方程式は

$$
\sum_{n=1}^{\infty}q(\alpha)^n=1.
$$

左辺は無限級数の和であるから，

$$
\frac{q(\alpha)}{1-q(\alpha)}=1.
$$

すなわち，

$$
q(\alpha)=\frac12.
$$

したがって，

$$
\left(\frac{\beta}{\beta+\alpha}\right)^k
=\frac12.
$$

より，マルサス係数は，

$$
\alpha=\beta\left(2^{1/k}-1\right).
$$

## シミュレーション

まずイメージを掴むために1回のシミュレーション結果を図示してみる．

![](/images/gamma_branching_process/tree.png)
*横軸がノードの発生時刻，縦軸が世代*

分枝過程がシミュレーションできていそうなことがわかる．

次に10000回のシミュレーション結果で，総ノード数の対数をプロットしてみる．

![](/images/gamma_branching_process/count.png)
*横軸がノードの発生時刻，縦ノード数の対数．オレンジの点線は傾き $\alpha$ の直線．*

傾きが上で求めた $\alpha$ と一致していそうなことがわかる．

R のコードはこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/gamma_branching_process.R