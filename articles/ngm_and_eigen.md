---
title: "次世代行列の固有値・固有ベクトルについて"
emoji: "🦔"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 線形代数]
published: false
---

## このノートについて

このノートは [西浦 博; 小林 鉄郎; 安齋 麻美; 合原 一幸; ナタリー・リントン. 感染症流行を読み解く数理（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) 第1章の「人口の異質性」の節の行間を補うためのものです．

p.8 に次のような命題が出ていますが証明は省略されています.  

1. 任意の初期感染者数（$I_{C,0}$ および $I_{A,0}$） に対し， 世代 $t$ の総感染者数に対する次世代 $t+1$ の総感染者数の比 $\lambda = (I_{C,t+1}+I_{A,t+1})/(I_{C,t}+I_{A,t})$ は一定の値に収束する．
2. 任意の初期感染者数（$I_{C,0}$ および $I_{A,0}$）に対し， 総感染者数に対する小児および成人の感染者数の比（$I_{C,t} /(I_{C,t}+I_{A,t})$ および $I_{A,t} /(I_{C,t}+I_{A,t})$）は一定の値に収束する．

ここでは私なりこの命題の数値例と証明を与えます．数値例は R によります．

## 数値計算してみる

まず命題中の記号の意味を説明する．$I_{C,t}$, $I_{A,t}$ はそれぞれ小児と成人の $t$ 世代めの感染者数を表す．毎回書くのが大変なので，$I_{C,t}$, $I_{A,t}$ をまとめて，

$$
I_t =
\begin{pmatrix}
I_{C,t}\\
I_{A,t}
\end{pmatrix}
$$

と置くことにする．また総感染者数は次のように $N_t$ と置く．

$$
N_t= I_{C,t}+I_{A,t}= \mathbf1^\top I_t,
$$

ここで $\mathbf1$ は要素がすべて 1 のベクトル

$$
\mathbf1=
\begin{pmatrix}
1\
1
\end{pmatrix}
$$

とした．

各成分が非負 $K_{ij} \ge 0$ の次世代行列 $K$ と $t$ 世代めの感染者数を用いて，$t+1$ 世代めの感染者数は次のようになるとする．

$$
I_{t+1} = K I_t.
$$

次世代行列 $K$ の $(i, j)$ 成分 $K_{ij}$ が一人あたりの感染者から生み出される2次感染者数を表している．

次のように次世代行列に具体的な数を与えて、$I_t$ を計算してみよう．

```r
K <- matrix(c(1.4, 0.4,
              0.6, 0.3), byrow = TRUE, nrow = 2)
```

行列を繰り返し掛ける関数を次のように宣言する．

```r
iter_nextstep <- function(n, K){
  I <- c(1, 0)
  I_hist <- matrix(0, nrow = 2, ncol = n+1)
  I_hist[, 1] <- I
  for(i in 1:n){
    I <- K%*%I
    I_hist[,i+1] <- I  
  }
  return(t(I_hist))  
}

n <- 10
I_hist <- iter_nextstep(n,K)
```

プロットしてみると次の図のように指数的増加する系列が得られたことがわかる．

![](/images/ngm_and_eigen/original1.png)
*世代ごとの感染者数の推移．Cが小児（Child），Aが成人（Adult）．*

命題にあるような比を取ってみよう．総感染者数 $N_t$ を次のように求める．

```r
sumI <- rowSums(I_hist)
```

$N_{t+1}/Nt$ すなわち `sumI[-1]/sumI[-length(sumI)]` をプロットしたものが次の図である．

![](/images/ngm_and_eigen/ratio11.png)
*世代ごとの2次感染者数の推移*

図から一定の値に収束していくことがわかる．この「一定の値」は実は次世代行列 $K$ の最大固有値である．

R では固有値・固有ベクトルは次のようにして求められる．

```r
> ei_K <- eigen(K)
> print(ei_K$values)
[1] 1.586546 0.113454
```

総感染者数のうちの小児・成人の感染者数の割合もプロットしてみよう．

![](/images/ngm_and_eigen/ratio21.png)
*小児・成人の感染者数に占める割合の推移*

図から，こちらも一定の値に収束していくことがわかる．この「一定の値」は次世代行列 $K$ の最大固有値に属する固有ベクトルを足して1になるよう正規化した値である．

```r
> norm_v <- ei_K$vectors[,1]/sum(ei_K$vectors[,1])
> print(norm_v)
[1] 0.6819585 0.3180415
```

もうひとつ別のパターンの次世代行列でもやってみよう．次のようにする．

```r
K2 <- matrix(c(0, 2,
              3, 0), byrow = TRUE, nrow = 2)
```

![](/images/ngm_and_eigen/original2.png)
*世代ごとの感染者数の推移．Cが小児（Child），Aが成人（Adult）．*

$N_{t+1}/Nt$ をプロットしたものが次の図である．

![](/images/ngm_and_eigen/ratio12.png)
*世代ごとの2次感染者数の推移*

図から，この場合は振動して一定の値に収束しないことがわかる．

総感染者数のうちの小児・成人の感染者数の割合もプロットしてみよう．

![](/images/ngm_and_eigen/ratio22.png)
*小児・成人の感染者数に占める割合の推移*

図から，やはりこの場合は振動して一定の値に収束しないことがわかる．

`K2` の固有値を見ると次のように最大固有値 $\rho$ が他の固有値 $\lambda$ の絶対値と等しくなっている（$\rho = \lvert\lambda\rvert$）．

```r
> ei_K2 <- eigen(K2)
> print(ei_K2$values)
[1]  2.44949 -2.44949
```

作図も含めた R のコード全体は以下に置く：
https://github.com/abikoushi/Zenn_content/blob/main/R/ngm_and_eigen.R

## 手計算してみる

### 準備

次世代行列 $K$ は各成分が非負 $K_{ij} \ge 0$ であった．加えて以下の仮定を置く．

1. 対角化可能
2. 最大固有値 $\rho$ は他の固有値 $\lambda$ の絶対値より大きい（$\rho > \lvert\lambda\rvert$）

もしかしたらもっと緩い条件で示せるかもしれないが，ここではこれでいく．

固有値を $\rho, \lambda$ とし，それぞれの固有値に属する固有ベクトルを $v,u$ とする．すなわち

$$
Kv=\lambda v,\quad Ku=\rho u, 
$$

を満たす．仮定 1 より固有ベクトル $u$, $v$ と，固有ベクトルを並べた行列の逆行列の列ベクトル $w$, $z$ を用いて

$$
\begin{pmatrix}w & z \end{pmatrix}K \begin{pmatrix}v & u \end{pmatrix}= 
\begin{pmatrix}
\rho & 0\\
0 & \lambda 
\end{pmatrix} 
$$

と対角化できる．$w' v=1$, $z' u=1$ かつ $w' u=0$, $z' v=0$ となるように正規化されている．

$v=(v_1, v_2)'$, $u=(u_1, u_2)'$ として成分を書き下すと，$K$ は次のように分解できる．

$$
\begin{aligned}
K &=
\begin{pmatrix}
v & u
\end{pmatrix}
\begin{pmatrix}
\rho & 0\\
0 & \lambda 
\end{pmatrix} 
\begin{pmatrix}
w & z
\end{pmatrix}\\
&= 
\begin{pmatrix}
v_1 & u_1\\
v_2 & u_2
\end{pmatrix}
\begin{pmatrix}
\rho & 0\\
0 & \lambda 
\end{pmatrix} 
\begin{pmatrix}
w_1 & z_1\\
w_2 & z_2\\
\end{pmatrix}\\
&= 
\begin{pmatrix}
\rho v_1 w_1 +\lambda u_1 w_2& \rho v_1 z_1 + \lambda u_2 z_1\\
\rho v_2 w_1 + \lambda u_2 w_2 & \rho v_2z_1+\lambda u_2z_2
\end{pmatrix}\\
& = \rho vw'+\lambda uz'.
\end{aligned}
$$

したがって

$$
K^t = \rho^tv w' + \lambda^tuz'.
$$

初期感染者数を $I_0$ とすると

$$
\begin{aligned}
I_t &= K^t I_0  \\
&=\rho^tv(w^\top I_0) + \lambda^tu(z^\top I_0).
\end{aligned}
$$

ここで $\alpha=w' x_0$, $\beta=z' x_0$ と置けば

$$
I_t=\alpha\rho^tv+\beta\lambda^tu.
$$

したがい，

$$
N_t=
\alpha\rho^t(\mathbf1^\top v)
+
\beta\lambda^t(\mathbf1^\top u).
$$

### 命題1

$N_{t+1}/N_t$ を求める．分子分母を $\rho^t$ で割ると，

$$
N_{t+1}/N_t = 
\frac{
\alpha\rho(\mathbf1^\top v)
+
\beta(\lambda/\rho)^t
\lambda(\mathbf1^\top u)
}
{
\alpha(\mathbf1^\top v)
+
\beta(\lambda/\rho)^t
(\mathbf1^\top u)
}.
$$

今，2つめの仮定より $|\lambda|<\rho$ なので $(\lambda/\rho) \rightarrow 0$ ($t \to \infty$) より，

$$
\lim_{t\to\infty}
\frac{N_{t+1}}{N_t}
=
\rho.
$$


### 命題2

小児の集団の感染者割合を考える．

$$
\frac{I_{C,t}}{N_t} =\frac{
\alpha v_1
+
\beta(\lambda/\rho)^tu_1
}
{
\alpha(\mathbf1^\top v)
+
\beta(\lambda/\rho)^t(\mathbf1^\top u)
}.
$$

極限を取れば

$$
\lim_{t\to\infty}
\frac{I_{C,t}}{N_t}
=
\frac{v_1}{\mathbf1^\top v}.
$$

同様に

$$
\lim_{t\to\infty}
\frac{I_{A,t}}{N_t}
=
\frac{v_2}{\mathbf1^\top v}.
$$

おしまい．