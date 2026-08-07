---
title: "離散時間の分岐過程と感染症疫学へのその応用"
emoji: "🪾"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 確率過程]
published: false
---

## はじめに

## シミュレーションしてみる

平均が $R_0$ となるよう，ここでの負の二項分布の確率関数は次のようにパラメータ化される．

$$
\Pr(X=x) = \frac{\Gamma(k+x)}{x!\Gamma(k)}\left(\frac{R_0}{R_0+k}\right)^x \left(1+\frac{R_0}{k}\right)^{-k}
$$

この分布の平均は $R_0$，分散は $R_0+R_0^2/k$ となる

R の `*nbinom` 関数を使う際は，次のように置くと上記の確率関数と一致する．

$$
\texttt{size} \leftarrow k, \quad \texttt{prob} \leftarrow k/(k+R_0)
$$

$$
\Pr(Y=y)=\begin{cases}
\frac{1}{(1+R_0/k)^k} & y=1\\
\frac{\prod_{j=0}^{y-2}(j/k+y)}{y!}\left(\frac{k}{R_0+k}\right)^{ky} \left(\frac{R_0k}{R_0+k}\right)^{y-1} & y\ge2
\end{cases}
$$

## 導出してみる