---
title: "時間変換による非定常ポアソン過程のシミュレーション：べき乗則と指数則"
emoji: "⏳️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 確率過程]
published: true
---

## はじめに

[近江崇宏・野村俊一 『点過程の時系列解析』（共立出版）](https://www.kyoritsu-pub.co.jp/book/b10003181.html) には次のような記述がある（7.2『ポアソン過程のシミュレーション』の最後の段落）．

> 定常ポアソン過程をシミュレーションした後に，適切な時間変換を行うことで非定常ポアソン過程を得ることも考えられる．しかしながら，一般的にこの時間変換は解析的に得ることができないため，数値的な方法に頼らざるを得ず，効率のよい方法とは言えない．

裏を返せば，変換が解析的に得られる特殊ケースならば時間変換による非定常ポアソン過程（non-homogeneous Poisson process）のシミュレーションは効率のよい方法になり得る可能性がある（もちろん論理的には「裏」は必ずしも真でないが）．ここではそれをやってみたい．

[鎌倉稔成　21世紀の統計科学　第3章「生存時間・再発事象分析：理論と応用」5.3.3 他のモデルと例題](https://www.cirje.e.u-tokyo.ac.jp/research/reports/R-15-2.pdf)（PDF直リンク）を見ると，比較的扱いやすい非定常ポアソン過程のモデルとして，強度関数をべき関数にした非定常ポアソン過程や強度関数を指数関数とした非定常ポアソン過程が紹介されている．このノートではそれぞれを「べき乗則ポアソン過程」，「指数則ポアソン過程」と呼ぶことにして，時間変換によるシミュレーションの方法を考える．


## べき乗則ポアソン過程のシミュレーション

次の累積強度関数を持つ非定常ポアソン過程をべき乗則ポアソン過程と呼ぶことにする（ワイブル過程と呼ばれることもあるが，ワイブル分布の再生過程と紛らわしいかと思い，"power-law poisson process" を直訳した）.

$$
H(t) =\left( \frac{t}{\alpha} \right)^\beta \quad (\alpha>0, \  \beta>0). 
$$

非定常ポアソン過程にしたがう確率変数 $t_1$ は分布関数 $F(t)=1-\exp(-H(t))$ を持つ分布に従う．

$z_1, z_2, \ldots, z_n$ を一様乱数として，逆関数法により，$t_1$ は，

$$
t_1=\alpha(- \log (z_1))^{1/\beta}
$$ 

で得られる．

$t_k$ ($k > 1$) は条件付き確率分布

$$
 P(T>t_k|T>t_{k-1})= \frac{\exp\left( -( \frac{t_k}{\alpha} )^\beta \right) } {\exp \left(- ( \frac{t_{k-1}}{\alpha})^\beta \right)}
$$

に従うから，逆関数法により，

$$
t_k=\alpha \left(- \log (z_k) + \left(\frac{t_{k-1}}{\alpha}\right)^\beta  \right)^{1/\beta}．
$$

で得られる．

すなわちべき乗則ポアソン過程のシミュレーションは，次式によりできる．

$$
\begin{aligned}
t_1 &=\alpha(- \log (z_1))^{1/\beta}\\
t_k &=\alpha \left(- \log (z_k) + \left(\frac{t_{k-1}}{\alpha}\right)^\beta  \right)^{1/\beta}．
\end{aligned}
$$

もう少し簡単な形になるように変形してみる．$t_2$ に $t_1$ を代入すると，

$$
t_2 =\alpha \left(- \log (z_2) + \log (z_1)  \right)^{1/\beta}．
$$

同様にして $t_k$ に $t_{k-1}$ を代入すると，

$$
t_k =\alpha \left(- \sum_{i=1}^k\log (z_i)  \right)^{1/\beta}．
$$

すなわち，R による実装は，一例として次のようにできる：

```r
Weibull_process <- function(n, alpha, beta) {
  z <- runif(n)
  alpha * cumsum(-log(z))^(1 / beta)
}
```

100回のシミュレーション結果から，イベントの累積発生回数を図示してみよう．

![](/images/nhpp_pow_exp/nhpp_pow.png)
*シミュレーション結果．x軸が時間，y軸がイベントの累積発生回数，青い点線が累積強度関数を表す.*

## 指数則ポアソン過程

べき乗則ポアソン過程のほかに比較的扱いやすい非定常ポアソン過程として，強度関数を次のようにしたものがある．

$$
h(t) = \exp(a+bt)
$$

特定の名前は見つからなかったが，指数的増加を表せる非定常ポアソン過程として，ここでは指数則ポアソン過程と呼ぶことにする．

累積強度関数は $h(u)$ を区間 $[0,t]$ で積分して，次のようになる．

$$
H(t) = (\exp(at+b) - \exp(b))/a
$$

非定常ポアソン過程にしたがう確率変数 $t_1$ は分布関数 $F(t)=1-\exp(-H(t))$ を持つ分布に従う．

$z_1, z_2, \ldots, z_n$ を一様乱数として，逆関数法により，$t_1$ は，

$$
t_1=(\log(\exp(b) - a\log (z_1)) -b)/a
$$ 

で得られる．

$t_k$ ($k > 1$) は条件付き確率分布

$$
P(T>t_k|T>t_{k-1})= \frac{\exp( -  (\exp(at_k+b) - \exp(b))/a} {\exp (- (\exp(at_{k-1}+b) - \exp(b))/a )}
$$

に従うから，逆関数法により，次式で得られる．

$$
t_k= \left(\log (- a \log (z_k) + \exp(at_{k-1}+b)) -b \right)/a．
$$

すなわち，指数則ポアソン過程のシミュレーションは次式によりできる．

$$
\begin{aligned}
t_1 &=(\log(- a\log (z_1) + \exp(b)) -b)/a\\
t_k &= \left(\log (- a \log (z_k) + \exp(at_{k-1}+b)) -b \right)/a, \quad (k\ge2)．
\end{aligned}
$$

これは，

$$
\begin{aligned}
w_0 &= \exp(b)\\
w_k &= - a \log (z_k) + w_{k-1}
\end{aligned}
$$

と置けば，次の単純な漸化式で表せる．

$$
t_k= \left(\log (w_k) -b \right)/a, \quad (k=1,\ldots,n)．
$$

注意として，上の式は $a < 0$ の場合も含めると対数の中身が負になってしまうことがある．$a<0$ のとき，累積強度の上限が次のように有限であるから，ある程度時間が経過するとそれ以上イベントが発生しなくなる．

$$
\lim_{t\to\infty} H(t) = - \exp(b)/a
$$

このことを「次のイベントの発生は無限に先」と捉えると R による実装は，一例として次のようになる：

```r
NHPP_exp <- function(n, a, b) {
  mlogz <- -log(runif(n))
  z <- exp(b) + a * cumsum(mlogz)
  t <- rep(Inf, n)
  idx <- z > 0
  t[idx] <- (log(z[idx]) - b) / a
  t
}
```

100回のシミュレーション結果から，イベントの累積発生回数を図示してみよう．

![](/images/nhpp_pow_exp/nhpp_exp.png)
*シミュレーション結果．x軸が時間，y軸がイベントの累積発生回数，青い点線が累積強度関数を表す.*

作図に用いた R のコードは以下にまとめて置く：

https://github.com/abikoushi/Zenn_content/blob/main/R/nhpp_powerlaw.R