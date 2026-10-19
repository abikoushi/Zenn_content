---
title: "logsumexp と softmax"
emoji: "🍦"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 数値計算]
published: true
---

## logsumexp

対数尤度の計算などで，$\log \sum \exp(\cdots)$ という，対数関数の中に足し算が入った形の式が出てくることがあります．このとき指数関数 $\exp$ の中に大きい数が入ると極端に大きい数となり，扱いにくいことがあります．そこで次のような変形を考えます．

$$
\begin{aligned}
\log\left[ \sum_{i=1}^n \exp(v_i) \right]
&=\log\left[ \sum_{i=1}^n \exp(m)\exp(v_i - m)\right] \\
&=\log\left[ \exp(m)\sum_{i=1}^n \exp(v_i - m)\right] \\
&= m + \log\left[\sum_{i=1}^n \exp(v_i - m)\right] .
\end{aligned}
$$

この関係を使って適当に大きい数 $m$ を引いておけば指数関数の中に大きい数が入るのを避けられます．$m$ はどうやって選べばいいんだという問題がありますが，$v_1,\ldots ,v_n$ の最大値としておけば確実です．


## softmax

多項分布を使ったモデル化などで，次のように総和が1になるようベクトルを正規化することがあります．

$$
s_i = \frac{\exp(v_i)}{\sum_{i=1}^n \exp(v_i)}.
$$

ベクトル $v=(v_1, \ldots, v_n)'$ に対してベクトル $s=(s_1, \ldots, s_n)'$ を対応させる関数をソフトマックス（softmax）関数と呼びます．

ソフトマックス関数についても指数関数の中に大きい数が入ることを避けるため，logsumexp のときと同様に，次のような変形を考えることができます．

$$
\begin{aligned}
s_i &= \frac{\exp(v_i)}{\sum_{i=1}^n \exp(v_i)} \\
&= \frac{\exp(m)\exp(v_i-m)}{\exp(m)\sum_{i=1}^n \exp(v_i-m)}\\
&= \frac{\exp(v_i-m)}{\sum_{i=1}^n \exp(v_i-m)}.
\end{aligned}
$$

## R による数値例

それぞれ次のように実装できます．

```r
logsumexp <- function(v){
  m <- max(v)
  m + log(sum(exp(v-m)))
}


softmax <- function(v){
  vm <- v - max(v)
  exp(vm)/sum(exp(vm))
}
```

まず小さい数を入れてみると，素朴な実装と同じ結果が得られることがわかります．

```r
> v1 <- c(1, 1, 2)
> 
> print( logsumexp(v1) )
[1] 2.551445
> print( log(sum(exp(v1))) ) 
[1] 2.551445
> 
> print( softmax(v1) )
[1] 0.2119416 0.2119416 0.5761169
> print( exp(v1)/sum(exp(v1)) ) 
[1] 0.2119416 0.2119416 0.5761169
```

大きい数を入れてみると，素朴な実装より数値的に安定していそうなことがわかります．

```r
> v2 <- c(10000, 10000, 10001)
> 
> print( logsumexp(v2) )
[1] 10001.55
> print( log(sum(exp(v2))) )
[1] Inf
> 
> print( softmax(v2) )
[1] 0.2119416 0.2119416 0.5761169
> print( exp(v2)/sum(exp(v2)) ) 
[1] NaN NaN NaN
```

おしまい．