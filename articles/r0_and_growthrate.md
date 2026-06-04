---
title: "基本再生産数と内的増殖率"
emoji: "📈"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学]
published: true
---

## このノートについて

このノートは [西浦博　編著『感染症流行を読み解く数理』（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) の第2章，第3章（ともに西浦博による）の一部を，自分なりに行間を埋めながら読んでみたもので [SIRモデルと基本再生産数](https://zenn.dev/abe2/articles/sir_model_and_r0) の続き的な内容です．

ただし，『感染症流行を読み解く数理』だけだとSIRモデルを元に定義された基本再生産数と一般的な基本再生産数（後述）の関係がよくわからなかったので，[Ma, J. (2020). Estimating epidemic exponential growth rate and basic reproduction number. Infectious Disease Modelling, 5, 129-141.](https://www.sciencedirect.com/science/article/pii/S2468042719300491) も参考にしました．

数値例は R 言語によります．


## 準備

単位時間あたりの変化量が自分自身に比例する関係を次のように微分方程式で表してみる．

$$
\frac{d}{dt} C(t) = rC(t). 
$$

この解は次の指数関数である．

$$
C(t) = C(0)\exp(rt). \tag{1}
$$

これが解になっていることは両辺を微分してみれば確かめられる． 今回は $r$ を内的増殖率と呼ぶことにする．

さて，感染症の感染者数を考えると，感染者数が多いほど感染する機会が多いため，流行初期は指数的増加に近い傾向になると想像できる．

## SIRモデルと内的増殖率の関係

SIRモデルは感染症の流行過程を記述するモデルで，次の連立微分方程式で表される．

$$
 \begin{aligned}
 {\frac {d}{dt}}S(t)&=-\beta S(t)I(t) \\
 {\frac {d}{dt}}I(t)&=\beta S(t)I(t)-\gamma I(t)\\
 {\frac {d}{dt}}R(t)&=\gamma I(t)
 \end{aligned}
$$

パラメータの置き方は [SIRモデルと基本再生産数](https://zenn.dev/abe2/articles/sir_model_and_r0) と同じなので説明は省略する．

いま，流行初期で $S(t)$ がほぼ1の状況を考えると，$I(t)$ についての方程式から次が成り立つ．

$$
 {\frac {d}{dt}}I(t) \approx (\beta -\gamma) I(t)
$$

この場合の解は，次の指数関数である．

$$
I(t) = I(0)\exp( (\beta -\gamma) t ).
$$

この式を (1) と比較すると，内的増殖率は $r = \beta-\gamma$ である．

SIRモデルの解とこの指数関数を並べてプロットしてみよう．

![](/images/r0_and_growthrate/SIR_expgrowth.png)
*黒い実線がSIRモデルのI，青い点線が指数関数． SIRモデルのI以外の解はグレーで表した．*

## 基本再生産数と内的増殖率の関係

2章では，感染後 $\tau$ 日目の感染者が生み出す2次感染率 $A(\tau)$ を用いて，次式で $R_0$ が定義がされる．

$$
R_0 = \int_{0}^{\infty} A(\tau) \, d\tau.\tag{2}
$$

その意味は，右辺より，1人あたりの感染者が生み出す2次感染者の全期間に渡る合計といったところである．$i(t)$ を時刻 $t$ における新規感染者数とし，次の積分方程式を考える．

$$
i(t) = \int_0^{\infty} i(t-\tau)A(\tau) \, d\tau.
$$

これは再生産方程式と呼ばれる．各時点の感染者それぞれが，それぞれ $A(\tau)$ だけ生み出した新たな感染者を合計している．流行初期で $i(t) = i(0)\exp( rt )$ が成り立っているとすると，

$$
i(0)\exp( rt ) = \int_0^{\infty}  i(0)\exp( r(t-\tau))A(\tau) \, d\tau, 
$$

より，

$$
1= \int_0^{\infty}  \exp( -r\tau)A(\tau) \, d\tau.
$$

$g(\tau) = A(\tau)/R_0$ と置くと，次を得る．

$$
1= R_0 \int_0^{\infty}  \exp( -r\tau)g(\tau) \, d\tau.
$$

$R_0$ と内的増殖率 $r$ の間に次の関係が成り立つ．

$$
R_0=\frac{1}{\int_0^{\infty}  \exp( -r\tau)g(\tau) \, d\tau}. \tag{3}
$$

右辺の分母の積分はちょうど $g(\tau)$ のモーメント母関数に $-r$ を代入したものである．また $g(\tau)$ は全区間に渡る積分が 1 になるように $A(\tau)$ を正規化したものであるから，確率密度関数である．$g(\tau)$ を世代時間の密度関数と呼ぶことがあるようだ．


### 指数分布の場合

いま，SIRモデルに戻って，感染者はレート $\gamma$ で除外されていくとする．さらに，よく使われる正の連続型の分布として一番かんたんそうな指数分布を考え，感染者はレート $\gamma$ の指数分布に従って除外されるとする．指数分布の分布関数を $F(\tau)$ と書くことにする．$A(\tau)$ は感染者が除外されずに残る確率に比例するとすると，

$$
A(\tau) \propto 1-F(\tau)= \exp(-\gamma \tau).
$$

$[0,\infty]$ での積分が 1 になるよう $A(\tau)$ を正規化しなおすと，次が得られる

$$
g(\tau) = \gamma\exp(-\gamma \tau).
$$

これはパラメータ$\gamma$ の指数分布の密度関数 $g(\tau)$ である．すなわち，$g(\tau)$ は再び指数分布である．レートパラメータ $\gamma$ の指数分布のモーメント母関数は，次式で表される．

$$
M(t) = \frac{\gamma}{\gamma-t}.
$$

SIRモデルの時の内的増殖率 $r=\beta-\gamma$ を用いると，(3)より，

$$
\begin{aligned}
R_0 &= \frac{1}{M(-r)} = \frac{1}{M(-(\beta-\gamma))}=\frac{1}{\gamma/\beta}\\
&= \beta/\gamma .
\end{aligned}
$$

これは『感染症流行を読み解く数理』の1章でSIRモデルに基づいてとして導入された $R_0=\beta/\gamma$ と一致する（[SIRモデルと基本再生産数](https://zenn.dev/abe2/articles/sir_model_and_r0) 参照）．(2) の方がより一般的な状況を考えていると言えるだろう．


### ガンマ分布の場合

一般の分布の場合，基本再生産数 $R_0$ は内的増殖率だけでなく，世代時間の分布 $g(\tau)$ に依存する．

ためしに $g(\tau)$ をガンマ分布として計算してみよう．次のようにした．

```r
R0hat_gamma <- function(m,CV){
  r <- 0.135
  ab2 <- (m*CV)^2
  b <- ab2/m
  a <- m/b
  int <- integrate(function(x){exp(-r*x)*dgamma(x,shape = a,scale = b)},0,Inf)
  return(1/int$value)  
}
```

ガンマ分布のモーメント母関数も本当は閉じた形で求まるが，ここでは `integrate` 関数の数値積分を使って評価している．

![](/images/r0_and_growthrate/fig3_3.png)
*基本再生産数のヒートマップ．横軸を世代時間の分布の平均，縦軸を変動係数とした．*

この図は『感染症流行を読み解く数理』の図3.3の再現である．世代時間の平均が長くて分散が小さいほど，基本再生産数は大きくなることがわかる．

作図に用いた R のコード全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/r0_and_growthrate.R