---
title: "グラフ上の離散化された拡散方程式"
emoji: "🔥"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [偏微分方程式, R, 数値計算]
published: true
---

## 準備

[離散化された拡散方程式の状態空間モデル](https://zenn.dev/abe2/articles/10bc59bec3280c) の続き的な内容．

1次元→2次元→3次元と拡張していってもいいのだけれど，ここではグラフ（ノードやエッジがある方のグラフ）上の拡散を考えてみる．

さて，1次元の拡散方程式グラフ上

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

![丸と矢印のグラフ](/images/diffeq_onthe_graph/diffeq1.jpg)

グラフの丸をノード，線をエッジと呼ぶことにする．

(1)式の右辺の $u_{t,x}$ の係数は $u_{t,x}$ の持つノードの数で決まっていることがわかる.

一方, $u_{t,x}$ の周辺のノード $u_{t,x+1}$ や $u_{t,x-1}$ の係数は隣接関係があれば 1，なければ 0 である．

変数を行列・ベクトルでまとめて次の図のように表せる．

![丸と矢印のグラフと行列](/images/diffeq_onthe_graph/diffeq2.jpg)

結果，時点 $t$ の状態（$U_t$）から次の $t+1$ での状態（$U_{t+1}$）を作る式

$$
U_{t+1}  =A U_t\\
$$

は，隣接関係を表す行列（Adjacency），エッジの数を表す行列（Degree）と単位行列（Unit）を使って，

$$
U_{t+1} = ( \text{Unit} +\{ \text{Adjacency} - \text{Degree}\})\cdot U_t
$$
 
のように表せることがわかった．

これをグラフ上の拡散の定義だと思って，隣接関係を次の図のような上下左右として数値計算してみよう．

![格子上の隣接関係](/images/diffeq_onthe_graph/diffeq3.jpg)


## Rによる数値例

といってもちょっと力尽きてしまって，いくつか図を描くだけである．

隣接関係を表す行列を作る関数は次のようにした．

```r
make_adjmat_2d <-function(nr, nc){
  nn = as.integer(nr*nc)
  A <- Matrix::spMatrix(nrow = nn, ncol = nn)
  edges = nr*(0:nc) 
  for(i in 1L:nn){
    if(i-nr > 0){ #left bound
      A[i,i-nr] <- 1
    }
    if( i %% nr != 1L ){ #top bound
      A[i,i-1] <- 1
    }
    if((i+nr) <= nn){ #right bound
      A[i,i+nr] <- 1
    }
    if( (i %% nr != 0L) & i+1L <= nn ){ #bottom bound
      A[i,i+1] <- 1
    }
  }
  return(A)
}
```

if 文を使って上下左右の端っこのところで隣接関係がないようにしている．

この行列をプロットするとこんなふうだ：

![](/images/diffeq_onthe_graph/adjmat.jpg)

1期先の状態を作るのは次のようにした．

```r
difit <- function(coef, u, nt){
  U <- matrix(0, nt+1, length(u))
  U[1,] <- u
  for(i in 1:nt){
    u2 <- coef%*%u
    U[i+1,] <- as.vector(u2)
    u <- u2
  }
  return(U)
}

```

これは `nt` ステップ文 for 文で回すだけであとは上の式そのままという感じ．

係数行列（上の式の $A$）は次のようにした．


```r
makecoef_exp <- function(r, A){
  diag(nrow = nrow(A)) + r*(A-diag(rowSums(A)))
}
```

これも上の式そのままという感じ．

この行列をプロットするとこんなふうだ：

![](/images/diffeq_onthe_graph/coefmat.jpg)

初期状態で2箇所だけ濃度が高いところを作って24ステップ回した．

```r
uini <- numeric(length = nrow(A))
uini[1] <- 1
uini[8] <- 1

U <- difit(beta, uini, 24)

print(rowSums(U))
```

次のように，全体での量が変化していないことをたしかめた．

```r
> print(rowSums(U))
 [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
```

拡散していく様子のGIFアニメ：

![](/images/diffeq_onthe_graph/tiles.gif)


いくつか切り口を変えた図も描いてみる．

とりあえず全部並べる：

![](/images/diffeq_onthe_graph/tiles.png)


横（x）方向に切り取って濃度の折れ線グラフ：

![](/images/diffeq_onthe_graph/xoriented.png)


縦（y）方向に切り取って濃度の折れ線グラフ：

![](/images/diffeq_onthe_graph/yoriented.png)


最後に使用したRのコード全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/28b54a4fa9641c360d8df0fb285b29e8d5b9f78c/R/diffusion_onthe_graph.R