---
title: "Diekmann et al. (2010) より，疫学における区画モデルと次世代行列の関係"
emoji: "💨"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 微分方程式, 線形代数]
published: TRUE
---

## このノートについて

下記の論文を読んで自分なりに理解したことのメモです．数理的に厳密な議論というよりは例を通じた理解を目指します．

Diekmann, O., Heesterbeek, J. A. P., & Roberts, M. G. (2009). The construction of next-generation matrices for compartmental epidemic models. Journal of the royal society interface, 7(47), 873. https://pmc.ncbi.nlm.nih.gov/articles/PMC2871801/

数値例は R によります．

## 例1：SIRモデル

次の連立微分方程式で表されるモデルを SIR モデルと呼ぶ．

$$
\begin{aligned}
\dot{S} &= -\beta \frac{SI}{N},\\
\dot{I} &= \beta \frac{SI}{N} - \gamma I,\\
\dot{R} &= \gamma I.
\end{aligned}
$$

SIRモデルには感受性者 $S$，感染期 $I$，回復・免疫保持 $R$ の3つの区画がある．各パラメータは次のような意味を持つ．

- $N$: 総人口
- $\beta$: 感染期の個体による感染伝播率
- $\gamma$: 回復率（$1/\gamma$ が平均感染期間）

パラメータと初期値を適当に選んで解を示す．

```r
library(deSolve)
pars  <- c(beta=0.3, gamma=0.1)
times <- seq(0, 100, by = 0.1)
ini = c(S=999, I=1, R=0)
ode_out <- ode(y=ini, times=times, func=SIRmod, parms=pars)
```

![](/images/compartment_ngm/SIR1.png)

$S \approx N$ が成り立つ感染流行初期を考える．$S = N$ が成り立つところを Diekmann et al. (2010) では無感染定常状態（infection-free steady state）と呼んでいる．感染流行初期は無感染定常状態の近傍である．このとき，感染 $I$ は次のようになる．

$$
\dot{I} = \beta I - \gamma I
$$

新規感染の発生に関する項のみに注目し，その係数を $T$ で表すことにすると

$$
T = \beta.
$$

これは1人の感染期個体 $I$ が単位時間あたりに生み出す新規感染者の平均を表している．

次に，新規感染以外の状態の変化（回復や死亡など）に関する項のみに注目し，その係数を $\Sigma$ で表すことにすると，

$$
\Sigma = -\gamma.
$$

初期の感染は時間 $t$ について指数関数 $\exp((T+\Sigma)t)$ で近似できる．この指数関数と先程の SIR モデルの解を重ねてみる．

![](/images/compartment_ngm/SIR_plusexp.png)
*点線が指数関数，実線が SIR モデルの解．*

この $T$ の逆数の符号反転 $-T^{-1} = 1/\gamma$ は，感染状態にある個体が（死亡または回復するまでに）感染状態で過ごす平均滞在時間という意味を持つ．

この分解に基づき $K$ を次のように定義する．

$$
K = -T \Sigma^{-1}
$$

$K= \beta/\gamma$ は「感染者が各状態で過ごす平均時間」と「その状態にいるときに新規感染を発生させる速度」を掛け合わせることで，「1人の感染者が生み出す新規感染者の総数」を表している．この量は基本再生算数と呼ばれる．

基本再生算数を $R_0=\beta/\gamma$ と置くと，$S(0) > N R_0$ のとき感染が指数関数的に広がり，$S(0) < N R_0$ のとき減衰する．

初期値を次のように振って SIR モデルの解より感染 $I$ についてプロットしてみよう．

```r
df_state <- data.frame(S=seq(1, 991, by=10)) %>% 
  mutate(I=10) %>% 
  mutate(R=1000-(S+I))
```

![](/images/compartment_ngm/SIR_thres.png)

初期値による色分けに着目すると，ある値を境に増加するか減少するかが分かれることがわかる．この「ある値」は基本再生算数により決まる．

## 例2：SEIRモデル

$$
\begin{aligned}
\dot{S} &= -\beta \frac{SI}{N},\\
\dot{E} &= \beta \frac{SI}{N} - \nu E,\\
\dot{I} &= \nu E - \gamma I,\\
\dot{R} &= \gamma I.
\end{aligned}
$$

SEIRモデルには感受性者 $S$，潜伏期 $E$，感染期 $I$，回復・免疫保持 $R$ の4つの区画がある．

各パラメータは次のような意味を持つ．

- $N$: 総人口
- $\beta$: 感染期の個体による感染伝播率
- $\nu$: 潜伏期から感染期への移行率（$1/\nu$ が平均潜伏期間）
- $\gamma$: 回復率（$1/\gamma$ が平均感染期間）

パラメータと初期値を適当に選んで解を示す．

![](/images/compartment_ngm/SEIR1.png)

感染している状態を表す変数は $E$ と $I$ である．Diekmann et al. (2010) では，このように区画モデルから感染している状態を表す変数のみを抜き出したものを感染サブシステム（infection subsystem）と呼んでいる．

無感染定常状態 $S=N$ の近傍で感染サブシステムは次のようになる．

$$
\begin{aligned}
\dot{E} &= \beta I - \nu E \\
\dot{I} &= \nu E - \gamma I
\end{aligned}
$$

$x = (E, I)'$ とまとめて置く．これを行列の和 $\dot{x} = (T + \Sigma)x$ の形に分解しよう．

SIR モデルのときと同様に，新規感染に関する項のみに注目した係数を $T$ とすると，

$$
T = \begin{pmatrix}
0 & \beta \\
0 & 0
\end{pmatrix}.
$$

このように感染の状態が複数ある一般の場合には $T$ は行列となる．Diekmann et al. (2010) では $T$ を伝播行列（transmission matrix）と呼ぶ．

感染以外の状態の変化に関する項のみに注目した係数を $\Sigma$ とすると，

$$
\Sigma = \begin{pmatrix}
-\nu  & 0 \\
\nu & -\gamma
\end{pmatrix}.
$$

Diekmann et al. (2010) では $\Sigma$ を推移行列（transision matrix）と呼ぶ．

初期の感染は時間 $t$ について行列の指数関数 $\exp((T+\Sigma)t)$ で近似できる．この指数関数と先程の SEIR モデルの解を重ねてみる．

![](/images/compartment_ngm/SEIR_plusexp.png)
*点線が行列の指数関数，実線が SEIR モデルの解．*

各状態で過ごす平均時間の行列 $-\Sigma^{-1}$ を求めると次のようになる．

$$
-\Sigma^{-1} = \begin{pmatrix}
1/\nu & 0 \\
1/\gamma & 1/\gamma
\end{pmatrix}
$$

この分解に基づき，$K_L$ を次で定義する．

$$
K_L = -T \Sigma^{-1}
$$

要素を計算すると次のようになる．

$$
K_L = \begin{pmatrix}
\beta /\gamma & \beta/\gamma \\
0 & 0
\end{pmatrix}
$$

$K_L$ の固有値は $\beta /\gamma, 0$ であり，最大固有値は $\beta /\gamma$ である．

基本再生算数を $R_0=\beta/\gamma$ と置くと，$S(0) > N R_0$ のとき感染が指数関数的に広がり，$S(0) < N R_0$ のとき減衰する．

初期値を次のように振って SEIR モデルの解より感染 $E+I$ についてプロットしてみよう．

```r
df_state <- data.frame(S=seq(1, 991, by=10)) %>% 
  mutate(E=0,I=10) %>% 
  mutate(R=1000-(S+E+I))
```

![](/images/compartment_ngm/SEIR_thres.png)

初期値による色分けに着目すると，ある値を境に増加するか減少するかが分かれることがわかる．この「ある値」は基本再生算数により決まる．

Diekmann et al. (2010) では $K_L$ を大きいドメインを持つ次世代行列（next-generation matrix with Large domain）と呼んでいる．

SEIRモデルでは，$T$ の2行目がすべて 0 であるため，感染直後の状態（states-at-infection）は $E$ のみである．したがって，不要な情報を削ぎ落とした通常の次世代行列（original next-generation matrix） $K$ は単にスカラーとなる．

$$
K = \frac{\beta}{\gamma}.
$$

通常の次世代行列 $K$ の場合も，（1行1列の行列と見れば）最大固有値は $\beta /\gamma$ である．

## まとめ

区画モデルの感染サブシステム $x$ を伝播行列 $T$ と 推移行列 $\Sigma$ を用いて $\dot{x} = (T + \Sigma)x$ の形に表せるとき，大きいドメインを持つ次世代行列は $K_L = -T \Sigma^{-1}$ で定義される．基本再生算数は $K_L$ の最大固有値で定義できる．

さらに，感染直後の状態のみに注目すると通常の次世代行列 $K$ を作れる．通常の次世代行列は大きいドメインを持つ次世代行列以下の次元を持つ．$K$ の最大固有値は $K_L$ の最大固有値と一致する．

## 関連記事など

駆け足になってしまった部分もあると思うので，わかりにくかった箇所は適宜以下の記事も参照してもらえるとうれしい．

- [SIRモデルと基本再生産数](https://zenn.dev/abe2/articles/sir_model_and_r0)
- [行列の指数関数と常微分方程式についてのイントロ](https://zenn.dev/abe2/articles/intro_matexp)

作図も含めたコード全体は以下に置く：

https://github.com/abikoushi/Zenn_content/blob/main/R/compartment_ngm.R