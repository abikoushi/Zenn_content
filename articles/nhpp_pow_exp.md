---
title: "時間変換による非定常ポアソン過程のシミュレーション：べき乗則と指数則"
emoji: "🌱"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R]
published: false
---

## はじめに

近江崇宏・野村俊一 『点過程の時系列解析』（共立出版）には次のような記述がある（7.2『ポアソン過程のシミュレーション』の最後の段落）

> 定常ポアソン過程をシミュレーションした後に，適切な時間変換を行うことで非定常ポアソン過程を得ることも考えられる．しかしながら，一般的にこの時間変換は解析的に得ることができないため，数値的な方法に頼らざるを得ず，効率のよい方法とは言えない．

裏を返せば，変換が解析的に得られる特殊ケースならば時間変換による非定常ポアソン過程のシミュレーションは効率のよい方法になり得る可能性がある（もちろん論理的には「裏」は必ずしも真でないが）．ここではそれをやってみる．

## ワイブル過程のシミュレーション

ワイブル過程（Weibull process）とは，累積強度関数 

$$
H(t) =\left( \frac{t}{\alpha} \right)^\beta
$$

を持つ非定常ポアソン過程（non homogeneous）の一種.

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

ワイブル過程のシミュレーションは，次の式でできた．

$$
\begin{aligned}
t_1 &=\alpha(- \log (z_1))^{1/\beta}\\
t_k &=\alpha \left(- \log (z_k) + \left(\frac{t_{k-1}}{\alpha}\right)^\beta  \right)^{1/\beta}．
\end{aligned}
$$

$t_2$ に $t_1$ を代入すると，

$$
t_2 =\alpha \left(- \log (z_2) + \log (z_1)  \right)^{1/\beta}．
$$

同様にして $t_k$ に $t_{k-1}$ を代入すると，

$$
t_k =\alpha \left(- \sum_{i=1}^k\log (z_i)  \right)^{1/\beta}．
$$

すなわち，R による実装は：

```r
Weibull_process <- function(n, alpha, beta) {
  z <- runif(n)
  alpha * cumsum(-log(z))^(1 / beta)
}
```


## 指数則ポアソン過程

ワイブル過程のほかに比較的扱いやすく，指数的増加・減少を表せる非定常ポアソン過程として，強度関数を次のようにしたものがある．

$$
h(t) = \exp(a+bt)
$$

例えば [鎌倉稔成　21世紀の統計科学　第3章「生存時間・再発事象分析：理論と応用」5.3.3 他のモデルと例題](https://www.cirje.e.u-tokyo.ac.jp/research/reports/R-15-2.pdf)（PDF直リンク）を参照．

累積強度関数は $h(t)$ を初期条件 0 で積分して，次のようになる．

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

上記のコードは $a < 0$ の場合も含めると `log()` の中が負になってしまうことがある．

$a<0$ のとき，累積強度の上限が次のように有限であるから，ある程度時間が経過するとそれ以上イベントが発生しなくなる．

$$
\lim_{t\to\infty} H(t) = - \exp(b)/a
$$

これは考えてみれば当たり前のことであった．なので次のように実装すると良さそうである．


```r
NHPP_exp = function(n, a, b){
  mlogz = -log(runif(n))
  t <- rep(Inf, n)
  t[1] <- (log(a*mlogz[1] + exp(b))-b)/a
  if(n >= 2){
    for(i in 2:n){
      s <- a*mlogz[i] + exp(a*t[i-1]+b) 
      if(s<0){
        break
      }
      t[i] =(log(s) - b)/a
    }
  }
  return(t)
}
```

![](https://static.zenn.studio/user-upload/579da094b926-20260508.png)



強度関数を $h(t) = \exp(a+bt)$ とした非定常ポアソン過程のシミュレーションは次の式でできた．

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

すなわち，R による実装は：

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