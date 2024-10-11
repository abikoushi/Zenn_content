---
title: "離散化されたグラフ上の拡散方程式"
emoji: "👌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [偏微分方程式,R,数値計算]
published: false
---

[離散化された拡散方程式の状態空間モデル](https://zenn.dev/abe2/articles/10bc59bec3280c) の続き的な投稿だが, 先の投稿を予備知識として要求しない.

さて，1次元の拡散方程式

$$
\frac{\partial}{\partial t} u(t,x) = a \frac{\partial}{\partial x^2} u(t,x)
$$

を離散化し, $u_{t, x}$ ($t=0,1, \ldots ,m$, $x=1, \ldots ,n$) について次の差分方程式を考える.

$$
\begin{aligned}
u_{t+1,x}-u_{t,x} &= a ( u_{t,x+1}-u_{t,x} - (u_{t,x}-u_{t,x-1})  )\\
&=a ( u_{t,x+1}-2u_{t,x}+u_{t,x-1}) 
\end{aligned}
$$

右辺の $u_{t,x}$ の係数は隣り合うノードの数で決まっていることがわかる.

![](/images/diffeq_onthe_graph/diffeq1.jpg)



