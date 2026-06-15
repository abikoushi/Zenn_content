---
title: "打ち切りのあるときの致命割合"
emoji: "✂️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学]
published: false
---

## はじめに

## シミュレーション

次のようなデータ生成過程を設定し，シミュレーションを行う．

$$
\begin{aligned}
t^c_i &\sim \mathrm{NHPP}(\lambda(t))\\
t^d_i &= t^c_i+x_i\\
(x_i \mid d_i=1) &\sim F(\cdot)\\
d_i &\sim \mathrm{Bernoulli}(p)
\end{aligned}
$$

$t^d_i <t$ を満たす $t^c_i$ の数を $C_t$ ，$t^d_i <t$ を満たす $t^d_i$ の数を $D_t$ とする．また $\lambda_t$ を次のようにおく．

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

