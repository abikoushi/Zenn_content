---
title: "シニングによる非定常ポアソン過程のシミュレーション"
emoji: "🪚"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 確率過程]
published: false
---

シニング（thinning）による非定常ポアソン過程のシミュレーションのコードを R で書いたので，簡単なメモとして残しておこう思う．

アルゴリズムの解説は例えば以下の文献にある．

- [近江崇宏・野村俊一『点過程の時系列解析』共立出版](https://www.kyoritsu-pub.co.jp/book/b10003181.html)
- [Jochen Voss, "An Introduction to Statistical Computing: A Simulation-based Approach." Wiley. ](https://www.wiley.com/en-no/shop/general-introductory-statistics/an-introduction-to-statistical-computing-a-simulation-based-approach-p-9781118357729)

アルゴリズムを愚直に書くと次のようになる．

```r
get_NHPP <-function(Tmax, lambda, ...,
                    lambda2=NULL, maxit=10000){
  if(is.null(lambda2)){
    opt <- optimise(lambda, lower = 0, upper = Tmax, maximum = TRUE,
                    ...)
    lambda2 <- opt$objective
  }
  Ts = c()
  Tast = 0
  for(i in 1:maxit){
    Tast = Tast + rexp(1, lambda2)
    if(Tast>Tmax){break}
    if(runif(1) < lambda(Tast, ...)/lambda2){
      Ts =c(Ts,Tast)
    }
  }
  return(Ts)
}
```

累積強度関数（cumulative intensity function）の逆関数が簡単な形で求まるとシニングのありがたみが薄れるので，強度関数（intensity function）は次のような周期 $P$ の関数にしてみる．

$$
\lambda(t) = a(\cos(2\pi t/P) + 1). 
$$

 1 を足しているのは単に強度関数を非負（0 以上）にしたいためである．

累積強度関数は，

$$
\Lambda(t) = a\left(\frac{P}{2\pi}\sin(2\pi t/P) + t\right) 
$$

である．

R では次のように書いた．

```r
lambda <- function(x, a, P){a*(cos(x*2*pi/P)+1)}
Lambda <- function(x, a, P){a*((P/(2*pi))*sin(x*2*pi/P)+x)}
```

次のように強度関数を与えてやると区間 $[0, 50]$ で非定常ポアソン過程のシミュレーションができる．

```r
dat <- get_NHPP(50, lambda, a=2, P=20)
```

とりあえず 100 回シミュレーションを回してプロットしてみる．

![](/images/nhpp_thinning/nhpp_cos.png)
*上のパネルはある時点までの累積イベント数，赤い点線が累積強度関数を表す．下のパネルは時間ごとの平均イベント数, 赤い点線が強度関数を表す．*

与えた強度関数の点過程がシミュレーションできていそうなことがわかる．

ところで，パラメータ $\theta$ を持つ強度関数の非定常ポアソン過程の尤度関数は，次のように書ける．

$$
L(\theta) = \prod_{i=1}^n \lambda(t_i; \theta) \exp(- \Lambda(t_i; \theta))
$$

 $\theta = (a, P)$ とまとめて置き，次のように対数尤度関数を書いて推定してみよう．

```r
loglik <- function(par, y, lambda, Lambda, maxT){
  a <- exp(par[1])
  P <- exp(par[2])
  loglambda <- log(lambda(y, a, P))
  sum(loglambda) - Lambda(maxT, a, P)
}
```

対数尤度関数を `optim` で最大化するには次のようにする．

```r
opt1 <-optim(c(0,3), loglik, y=dat,
             lambda = lambda, Lambda = Lambda,
             maxT = 50,
             control = list(fnscale=-1),
             method = "Nelder-Mead")
```

次の図のように最尤推定されたパラメータを用いてプロットしてみると，手元のデータへの当てはまりがわかる．

![](/images/nhpp_thinning/nhpp_fit.png)

作図も含めたコード全体は以下に置く：

https://github.com/abikoushi/Zenn_content/blob/main/R/nhpp_thinning.R


関連記事：[時間変換による非定常ポアソン過程のシミュレーション：べき乗則と指数則](https://zenn.dev/abe2/articles/nhpp_pow_exp)