---
title: "ポアソン分布と多項分布およびガンマ分布とディリクレ分布の関係"
emoji: "🕌"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学, 確率分布]
published: true
---

## 概略

独立な $k$ 個のポアソン分布を和で条件付けると多項分布になる．

また，独立な $k$ 個のガンマ分布を和で条件付けるとディリクレ分布になる．


![](https://storage.googleapis.com/zenn-user-upload/f13d51385bf0-20260107.jpg)
*示したいことの概念図*

https://gist.github.com/abikoushi/46c4fc1a1185ec7b985c1692617fdd37

## ポアソン分布から多項分布へ

$X_j$ をポアソン分布に従う確率変数，

$$
X_j \sim \mathrm{Poisson}(\lambda_j)
$$

とする. $X_j$ ($j=1,\ldots, k$) の総和を $n$ として条件付ける．

$$
\begin{aligned}
&P\!\left(X_j = x_j \;\middle|\; \sum_{i=1}^k X_i = n \right)\\
&= \frac{ P\!\left(X_j = x_j \right) \; P\!\left( \sum_{i=1}^k X_i = n \;\middle|\; X_j = x_j \right) }
       { P\!\left( \sum_{i=1}^k X_i = n \right) }\\
&= \frac{ P\!\left(X_j = x_j \right) \; P\!\left( \sum_{i \ne j} X_i = n - x_j \right) }
       { P\!\left( \sum_{i=1}^k X_i = n \right) }\\
&= \frac{\mathrm{Poisson}(x_j|\lambda_j) \cdot \mathrm{Poisson}(n-x_j|\sum_{i\neq j} \lambda_j)}{\mathrm{Poisson}(n|\sum_{i=1}^k \lambda_i)}\\
&= \frac{ \dfrac{\lambda_j^{x_j}}{x_j!} e^{-\lambda_j} 
        \cdot \dfrac{ \left( \sum_{i \ne j} \lambda_i \right)^{\,n-x_j} }{(n-x_j)!} 
        e^{ - \sum_{i \ne j} \lambda_i } }
       { \dfrac{ \left( \sum_{i=1}^k \lambda_i \right)^n }{n!} 
        e^{ - \sum_{i=1}^k \lambda_i } }\\
& = \frac{n!}{x_j! (n-x_j)!} 
   \left( \frac{\lambda_j}{\sum_{i=1}^k \lambda_i} \right)^{x_j}
   \left( \frac{\sum_{i \ne j} \lambda_i}{\sum_{i=1}^k \lambda_i} \right)^{n-x_j}.
\end{aligned}
$$

これは二項分布の確率関数である．

すなわち，$\bar{\lambda}_j=\frac{\lambda_j}{\sum_j \lambda_j}$ とおくと，

$$
P\!\left(X_j = x_j \;\middle|\; \sum_{i=1}^k X_i = n \right) = \mathrm{Binomial}(n,\bar{\lambda}_j).
$$

これはすべての $j \;(j=1,\dots,k)$ で成り立つが，ここでは $x_1, x_2, \ldots, x_k$ と順に考えていくことにする．

$x_1, x_2, \ldots, x_k$ の同時分布は

$$
\begin{aligned}
&P\left(X_1=x_1,X_2=x_2, \ldots, X_k=x_k \;\middle|\; \sum_{i=1}^k=n\right)\\
&=P\left(X_k=x_k \;\middle|\; X_1=x_1, X_2 = x_2, \ldots, X_{k-1}=x_{k-1}, \sum_{i=1}^k=n\right)\\
&\cdots P\left(X_2=x_2 \;\middle|\; X_1=x_1, \sum_{i=1}^k=n\right)P\left(X_1=x_1 \;\middle|\; \sum_{i=1}^k=n\right)
\end{aligned}
$$

であるから，

$$
\begin{aligned}
&P\!\left(X_2 = x_2 \;\middle|\; \sum_{i=1}^k X_i= n, X_1 = x_1 \right)\\
&=P\!\left(X_2 = x_2 \;\middle|\; \sum_{i=2}^k X_i = n-x_1 \right),\\
&P\!\left(X_3 = x_3 \;\middle|\; \sum_{i=1}^k X_i= n, X_1 = x_1 \right)\\
&=P\!\left(X_3 = x_3 \;\middle|\; \sum_{i=3}^k X_i = n-(x_1+x_2) \right),\\
& ~ \vdots\\
&P\!\left(X_k = x_k \;\middle|\; \sum_{i=1}^k X_i= n, X_1 = x_1, \ldots, X_{k-1} = x_{k-1} \right)\\
&=P\!\left(X_k = x_k \;\middle|\; X_k = n-\sum_{i=1}^{k-1}x_i \right),
\end{aligned}
$$

より，

$$
\begin{aligned}
&P \!\left( X_1 = x_1, \ldots,  X_k = x_k \;\middle|\;\sum_{i=1}^k X_i = n \right)\\
 &=\frac{n!}{x_1! \cdots x_k!} \prod_{j=1}^k
   \left( \frac{\lambda_j}{\sum_{i=1}^k \lambda_i} \right)^{x_j}.
\end{aligned}
$$

これは多項分布の確率関数である．

## ガンマ分布からディリクレ分布

$X_j$ をガンマ分布に従う確率変数，

$$
X_j \sim \mathrm{Gamma}(\alpha_j, \beta)
$$

とする. $X_j$ ($j=1,\ldots, k$) の総和を 1 として条件付ける．

$$
\begin{aligned}
& \quad P\!\left(X_j = x_j \;\middle|\; \sum_{i=1}^k X_i = 1 \right)\\
& = \frac{ P\!\left(X_j = x_j \right) \; P\!\left( \sum_{i=1}^k X_i = 1 \;\middle|\; X_j = x_j \right) }
       { P\!\left( \sum_{i=1}^k X_i = 1 \right) } \\
& = \frac{ P\!\left(X_j = x_j \right) \; P\!\left( \sum_{i \ne j} X_i = 1 - x_j \right) }
       { P\!\left( \sum_{i=1}^k X_i = 1 \right) }\\
&= \frac{\mathrm{Gamma}(x_j|\alpha_j, \beta) \cdot \mathrm{Gamma}(1-x_j|\sum_{_{i\neq j}} \alpha_i, \beta)}{\mathrm{Gamma}(1|\sum_{i=1}^k  \alpha_i, \beta)}\\
& = \frac{ \dfrac{\beta^{\alpha_j}}{\Gamma(\alpha_j)} x_j^{\alpha_j - 1} e^{-\beta x_j}
        \cdot \dfrac{\beta^{\sum_{i \ne j} \alpha_i}}{\Gamma\!\left(\sum_{i \ne j} \alpha_i\right)}
         (1-x_j)^{\sum_{i \ne j} \alpha_i - 1} e^{ -\beta (1-x_j)} }
       { \dfrac{\beta^{\sum_{i=1}^k \alpha_i}}{\Gamma\!\left(\sum_{i=1}^k \alpha_i\right)} }\\
&= \frac{\Gamma\!\left(\sum_{i=1}^k \alpha_i\right)}{\Gamma(\alpha_j)\,\Gamma\!\left(\sum_{i \ne j} \alpha_i \right)}
   x_j^{\alpha_j - 1} (1-x_j)^{\sum_{i \ne j} \alpha_i - 1}
\end{aligned}
$$

これはベータ分布の密度関数である．すべての $j \;(j=1,\dots,k)$ で成り立つ．

多項分布のときと同様にして，同時分布は

$$
\begin{aligned}
&P \!\left( X_1 = x_1, \ldots,  X_k = x_k \;\middle|\;\sum_{i=1}^k X_i = 1 \right)\\
& = \frac{\Gamma\!\left(\sum_{i=1}^k\alpha_i\right)}{\Gamma(\alpha_1)\cdots\Gamma(\alpha_k) }
\prod_{j=1}^k x_j^{\alpha_j - 1} 
\end{aligned}
$$

これはディリクレ分布の密度関数である．

## 数値的確認

僕自身はここで述べたような確率分布に関する話の大半を「乱数づくりゲーム」のアナロジーで理解している．私見だがそういう人は結構多いのではないかと思う．

まずポアソン分布に従う（擬似）乱数をつかって多項分布に従う乱数を作ってみよう．

R では関数 `rpois(n , rate)`  でポアソン分布に従う乱数を作れる．

$n$ と $\lambda_j$ に適当な数を与えて，発生させた乱数から総和が $n$ になるものだけをとってくれば「総和を $n$ として条件付ける」に相当する．今回は $n=10$ としよう．

```r
poi_to_mult <- function(N,lambda){
  X = matrix(0L, N, length(lambda))
  i = 0L
  while (i<=N) {
    x = rpois(length(lambda),lambda)
    if(sum(x) == 10L){
      X[i,] = x
      i = i+1L
    }
  }  
  return(X)
}

lambda = c(1,2,4,3)
system.time({
  set.seed(1234);res_mult <- poi_to_mult(10000,lambda)  
})
#  user  system elapsed 
# 0.066   0.000   0.066 
```

多変量の分布の可視化はいつも悩むが今回は周辺分布だけプロットしてみることにする．ほかにもアイデアがあれば教えてほしい．

![](/images/poi-mul_gam-dir/mult_cdf.png)
*実践が経験分布．点線が理論的な（？）分布関数*

乱数の経験分布が上で求めた周辺分布の二項分布とほぼ一致していることがわかる．

ガンマ分布は連続型なので総和がぴったり1になる確率はゼロになってしまい，同様にやるのは難しい．

しかし，ガンマ分布に従う確率変数はスケール変換（定数倍）したものもガンマ分布に従う性質がある．

$$
X \sim \mathrm{Gamma}(a,1)
$$

ならば

$$
(bX) \sim \mathrm{Gamma}(a,b)
$$

である. すなわち，総和が1になるように独立なガンマ乱数をスケール変換してやるとディリクレ分布に従う乱数が得られる．

こんなふうにした：

```r
gamm_to_dir <- function(N, alpha){
  X = matrix(0L, N, length(alpha))
  for(i in 1:10000){
    x = rgamma(length(alpha), alpha, 1)
    X[i,] = x/sum(x)    
  }
  return(X)
}

alpha = c(1,2,4)
system.time({
  set.seed(2345); res_dir <- gamm_to_dir(10000, alpha)
})
#  user  system elapsed 
# 0.009   0.000   0.009 
```

![](/images/poi-mul_gam-dir/dir_cdf.png)
*実践が経験分布．点線が理論的な（？）分布関数*

乱数の経験分布が上で求めた周辺分布のベータ分布とほぼ一致していることがわかる．

作図も含めたコードの全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/poi-mul_gam-dir.R

おしまい．
