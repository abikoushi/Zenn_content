---
title: "離散化されたグラフ上の拡散方程式"
emoji: "👌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [偏微分方程式,R,数値計算]
published: false
---

[離散化された拡散方程式の状態空間モデル](https://zenn.dev/abe2/articles/10bc59bec3280c) の続き的な内容だが, 先の投稿を予備知識として要求しない.

さて，1次元の拡散方程式

$$
\frac{\partial}{\partial t} u(t,x) = a \frac{\partial}{\partial x^2} u(t,x)
$$

を離散化し, $u_{t, x}$ ($t=0,1, \ldots ,m$, $x=1, \ldots ,n$) について次の差分方程式を考える.

$$
\begin{aligned}
u_{t+1,x}-u_{t,x} &= a ( u_{t,x+1}-u_{t,x} - (u_{t,x}-u_{t,x-1})  )\\
&=a ( u_{t,x+1}-2u_{t,x}+u_{t,x-1}) \tag{1}
\end{aligned}
$$

この差分方程式を次の画像のようにグラフで描いて考えてみる．

![](/images/diffeq_onthe_graph/diffeq1.jpg)

グラフの丸をノード，線をエッジと呼ぶことにする．

(1)式の右辺の $u_{t,x}$ の係数は $u_{t,x}$ の持つノードの数で決まっていることがわかる.

一方, $u_{t,x}$ の周辺のノード $u_{t,x+1}$ や $u_{t,x-1}$ 





