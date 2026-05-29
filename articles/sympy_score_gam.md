---
title: "SymPyの記号計算を使ってスコア検定を作る：ガンマ分布の形状パラメータの検定を例に"
emoji: "📖"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, SymPy, 統計学]
published: true
---

## はじめに

この記事は [Sibley A, Li Z, Jiang Y, Li YJ, Chan C, Allen A, Owzar K. Facilitating the Calculation of the Efficient Score Using Symbolic Computing. Am Stat. 2018;72(2):199-205.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6092959/) にインスパイアされた．

ガンマ分布の密度関数は次式で与えられる．

$$
f(x;\alpha,\theta)= \frac{1}{\Gamma(\alpha)\theta^\alpha} x^{\alpha-1}e^{-x/\theta}, \qquad x>0.
$$

$\alpha$ を形状（shape）パラメータ， $\theta$ を尺度（scale）パラメータと呼ぶ．

ガンマ分布は形状パラメータによって大きく性質が変わる. 形状パラメータが 1 より大きいときはモードの周りに集中した分布になり，1 以下のときは 0 で発散する裾の重い分布になる.

![](/images/sympy_score_gam/density.png)
*ガンマ分布の密度関数*

仮にお金の分布だとすると形状パラメータが 1 以下のときは貧富の差が大きい社会と言えそうだ．

一方で尺度パラメータが $k$ 倍になることはガンマ分布に従う確率変数を $k$ 倍することと等しい．尺度パラメータは文字通りスケールの取り方に依存して変わりうる．

そこで尺度パラメータを局外パラメータ（nuisance parameter）として形状パラメータの検定を行うことを考えてみたのが [ガンマ分布の形状パラメータの検定（あるいは局外パラメータがあるときのスコア検定の例題）](https://zenn.dev/abe2/articles/gam_shp_score) である．

すなわち，$X_1,\dots,X_n$ を独立に形状パラメータ $\alpha>0$，尺度パラメータ $\theta>0$ のガンマ分布に従う確率変数とし，このときに，

- 興味のあるパラメータ：形状パラメータ $\alpha$
- 局外パラメータ：尺度パラメータ $\theta$

として，帰無仮説 $\alpha=\alpha_0$ を検定するスコア検定を導いた．

さて，このくらいの手計算ならなんとかなるが，もうちょっと複雑なモデルになると解析的に求まるとしてもミスが増えそうだ．そこで今回は，より複雑なモデルへの応用をにらんで，同じことを [SymPy](https://www.sympy.org/en/index.html) の記号計算を使ってやってみる．ただし私が R になれているという理由から，[reticulate](https://rstudio.github.io/reticulate/) で R を経由して SymPy を使う．

## 本題

まずは SymPy をインストールする．

```r
library(reticulate)
virtualenv_install("r-reticulate", "sympy")
```

すでにインストールされているならば次のようにインポートすれば使用できる．

```r
library(reticulate)
reticulate::py_require("sympy")
sympy <- import("sympy")
stats <- import("sympy.stats")
```

記号を宣言する．特に `X` は後で期待値を計算したいので次のようにする．

```r
alpha <- sympy$Symbol("alpha", positive = TRUE)
theta <- sympy$Symbol("theta", positive = TRUE)
X <- stats$Gamma("X", alpha, theta)
```

サンプルひとつあたりの尤度，すなわち密度関数の対数を宣言し，スコアとフィッシャー情報量を求める．

```r
logpdf <- sympy$simplify(sympy$log(stats$density(X)(X)))

# score
S_alpha <- sympy$simplify(sympy$diff(logpdf, alpha))
S_theta <- sympy$simplify(sympy$diff(logpdf, theta))

# Fisher information
I_11 <- sympy$simplify(-stats$Expectation(sympy$diff(S_alpha, alpha)))
I_12 <- sympy$simplify(-stats$Expectation(sympy$diff(S_alpha, theta)))
I_21 <- sympy$simplify(-stats$Expectation(sympy$diff(S_theta, alpha)))
I_22 <- sympy$simplify(-stats$Expectation(sympy$diff(S_theta, theta)))
```

有効情報量を求める．

```r
I_eff <- sympy$simplify((I_11 - I_12 * I_21 / I_22))
```

プリントすると [ガンマ分布の形状パラメータの検定（あるいは局外パラメータがあるときのスコア検定の例題）](https://zenn.dev/abe2/articles/gam_shp_score) で求めたものと一致していることがわかるはずだ．ただしここでは $n=1$ である．

```r
> print(I_eff)
polygamma(1, alpha) - 1/alpha
```

データは乱数で与えることにする．帰無仮説は $\alpha=1$ とする．

```r
# データ
set.seed(1234)
x_data <- rgamma(5, shape = 2)

#帰無仮説のalpha
alpha0 <- 1
```

スコアと有効情報量をそれぞれサンプルの数だけ足し算する．

```r
# スコアの和
S_alpha0 <- sympy$Integer(0)
S_theta0 <- sympy$Integer(0)
I_eff0  <- sympy$Integer(0)
for(val in x_data) {
  S_alpha0 <- S_alpha0 + S_alpha$subs(dict(alpha = alpha0, X = val))
  S_theta0 <- S_theta0 + S_theta$subs(dict(alpha = alpha0, X = val))
  I_eff0 <- I_eff0 + I_eff$subs(dict(alpha=alpha0))
}
```

$\theta$ についてのスコアを 0 とおいて解く．

```r
sol_theta <- sympy$solve(sympy$Eq(S_theta0, 0), theta)
```

$\hat \theta$ を有効スコアに代入して検定統計量を求める．

```r
U_eff0 <- sympy$simplify(
  S_alpha0$subs(dict(theta=sol_theta[[1]]))^2/I_eff0$subs(dict(theta=sol_theta[[1]]))
)
```

R の数値に変換する．

```r
builtins <- import_builtins()
U_eff0_val <- py_to_r(builtins$float(U_eff0$evalf()))
```

[ガンマ分布の形状パラメータの検定（あるいは局外パラメータがあるときのスコア検定の例題）](https://zenn.dev/abe2/articles/gam_shp_score) で求めた次の関数と比較してみよう．

```r
gamma_shape_score_test <- function(x, alpha0) {
  n <- length(x)
  xbar <- mean(x)
  logxbar <- mean(log(x))
  
  U <- n * (logxbar - log(xbar) + log(alpha0) - digamma(alpha0))
  Ieff <- n * (trigamma(alpha0) - 1 / alpha0)
  stat <- (U^2) / Ieff
  pval <- pchisq(stat, df = 1, lower.tail = FALSE)
  
  list(
    statistic = stat,
    p.value = pval,
    null.value = alpha0
  )
}

res_r <- gamma_shape_score_test(x = x_data, 1)
```

ほぼ同じ値が得られたことがわかるはずだ．

```r
> print(U_eff0_val)
[1] 1.145261
> print(res_r$statistic)
[1] 1.145261
> print(all.equal(res_r$statistic, U_eff0_val))
[1] TRUE
```

R のコード全体はこちら：
https://github.com/abikoushi/Zenn_content/blob/main/R/sympy_score_gam.R

おしまい．