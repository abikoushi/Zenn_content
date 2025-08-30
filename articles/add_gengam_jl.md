---
title: "Distributions.jl を利用して確率分布を宣言する（一般化ガンマ分布の例）"
emoji: "⛑️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [Julia, R]
published: true
---

## 主にやったこと

次の分布関数 $f(x;a,b,k)$ を持つ確率分布を一般化ガンマ分布（[Generalized gamma distribution - wikipedia](https://en.wikipedia.org/wiki/Generalized_gamma_distribution)）という.

$$
F(x; a, b, k) = \frac{\gamma(a/k, (x/b)^k)}{\Gamma(a/k)}, \quad x > 0. \tag{1}
$$

ここで $\gamma(a,x)$ は第1種不完全ガンマ関数（ [不完全ガンマ関数 - wikipedia](https://ja.wikipedia.org/wiki/不完全ガンマ関数) ），パラメータ $(a, b, k)$ はすべて正の値を取る．密度関数は次の通り．

$$
f(x; a, b, k) = \frac{(k/b^a) x^{a-1} \exp(-(x/b)^k)}{\Gamma(a/k)}, \quad x > 0.
$$

`Distributions.jl` はいい感じのパッケージだが，幸か不幸か，一般化ガンマ分布は用意されていないようなので （ [Distributions.jl - Univariate Continuous Distributions](https://juliastats.org/Distributions.jl/stable/univariate/#Continuous-Distributions) ）これを実装してみた．具体的には，以下のコードのようにした．

https://github.com/abikoushi/Zenn_content/blob/main/jl/add_gengam_jl.jl

- 1行目から95行目：追加したい機能を `module` ... `end` で囲んでいます．こうしておくと他の関数と衝突したりしないし，「やっぱりこの書き方やめた！」となったとき上書きできるので試行錯誤段階からこうしておくほうが便利です．
- 25行目：`gamma_inc` は全区間での積分が1になるよう正規化された不完全ガンマ関数……つまり(1)式の分母の部分も含めて定義されています．戻り値の1つめが第1種不完全ガンマ関数（下側），2つめが第2種不完全ガンマ関数（上側）です．
- 75行目から78行目：乱数生成 `rand` は逆関数法によります．
- 88行目から91行目：対数の平均（≠平均の対数）を求める関数を「幾何平均の対数」というニュアンスを込めて `loggeommean` という名前で宣言しています．
- 96行目以降：正しく実装できているか答えあわせのためにやっています．下の解説に続きます．

`RCall.jl` を使って R を呼び，R の実装と比較して答えあわせとする．ただし R の `ggamma` パッケージではパラメトライズが (1) と少し異なり，こんなふうだ．

$$
F(x; a, b, k) = \frac{\gamma(k, (x/a)^k)}{\Gamma(k)}, \quad x > 0.
$$

下記の箇条書きのように置き換えると一致する．

- $a \leftarrow b$
- $b \leftarrow k$
- $k \leftarrow k/a$

すなわち，こんなふうにした．

```r
dggamma2 = function(x,a,b,k){ggamma::dggamma(x,b,k,a/k)}
```

パラメータ $(a, b, k)$ の値をいくつか選んで密度関数をプロットしたものが下の図だ.

![](/images/add_gengam_jl/dggamma.png)

R の実装（なめらかな曲線）と Julia の実装（点）で値が一致していることがわかる．

疑似乱数の経験分布と真の分布関数を比較したのが下の図．

![](/images/add_gengam_jl/ecdf_ggamma.png)

平均（`mean`）と対数の平均（`loggeommean`）については次のように標本平均と近い値が得られていることを見ている．

```julia
mean(x) #6.012540422658787
GG.mean(ggdist1) #6.0
mean(log, x) #0.848834299052846
GG.loggeommean(ggdist1) #0.8455686701969344
```

## 付録：一般化ガンマ分布の期待値，対数の期待値

確率変数 $X$ が一般化ガンマ分布に従うことを $X \sim \mathrm{GG}(a, b, k)$ と書くことにする. $X \sim \mathrm{GG}(a, 1, k)$ の場合の期待値 $E[X]$ を定義通り計算する．

準備としてガンマ関数の定義を確認する．

$$
\Gamma(z) = \int^{\infty}_{0} t^{z-1}\exp(-t) \, dt
$$

である. 期待値は次のようになる．

$$
\begin{aligned}
E[X] &= \int^{\infty}_{0} x \frac{k x^{a-1} \exp(-x^k)}{\Gamma(a/k)} \, dx \\
 &= \frac{k}{\Gamma(a/k)} \int^{\infty}_{0} x^{a} \exp(-x^k)\, dx.
\end{aligned}
$$

ここで $y=x^k$ と置くと, $x=y^{1/k}$, $dx/dy = (1/k)y^{1/k-1}$ より, $x^{a}\,dx=(1/k)y^{(a+1)/k-1} \, dy$ なので,

$$
E[X] = \frac{\Gamma((a+1)/k)}{\Gamma(a/k)}
$$

であるから, $X \sim \mathrm{GG}(a, b, k)$ のときは 

$$
E[X]=b\frac{\Gamma((a+1)/k)}{\Gamma(a/k)}
$$

次に, $X$ の自然対数をとった確率変数 $\log X$ の期待値 $E[\log X]$ を求めたい（そういうのを求めたい状況がたまにある）. $u = \log(t)$ と変数変換すると， $du = dt/t$ より, ガンマ関数は，次のように書ける.

$$
\Gamma(z) = \int^{\infty}_{-\infty} \exp(uz-\exp(u)) \, du
$$

さて，一般化ガンマ分布のパラメータ$b$はスケールパラメータで，$X \sim \mathrm{GG}(a, 1, k)$ のとき， $X$ を $b$ 倍した $bX$ は $\mathrm{GG}(a, b, k)$ に従う．そこでまず $X \sim \mathrm{GG}(a, 1, k)$ の場合を考える.

$$
\begin{aligned}
E[\log X] &= \int^{\infty}_{0} \log (x) \frac{k x^{a-1} \exp(-x^k)}{\Gamma(a/k)} \, dx \\
 &= \frac{1}{\Gamma(a/k)} \int^{\infty}_{0} k\log (x) x^a \exp(-x^k)\, dx/x
\end{aligned}
$$

$y=k\log (x)$と置き, $dy/k=dx/x$より,

$$
\begin{aligned}
E[\log X] &= \frac{1}{\Gamma(a/k)} \int^{\infty}_{-\infty} y \exp((a/k)y-\exp(y))\, (dy/k)\\
&= \frac{1}{\Gamma(a/k)} \int^{\infty}_{-\infty} \frac{\partial}{\partial a} \exp((a/k)y-\exp(y))\, dy\\
&= \frac{1}{\Gamma(a/k)} \frac{\partial}{\partial a}\left( \int^{\infty}_{-\infty} \exp((a/k)y-\exp(y))\, dy\right)\\
&=\left(\frac{\partial}{\partial a}\Gamma(a/k)\right)\frac{1}{\Gamma(a/k)}\\
&=\psi(a/k)/k.
\end{aligned}
$$

ここで$\psi(z)$はディガンマ関数 


$$
\psi(z) = \frac{d}{dz} \log \Gamma(z)
$$

である.

$E[\log(bX)]=E[\log(X)]+\log(b)$ なので， $X \sim \mathrm{GG}(a, b, k)$ のときは 

$$
E[\log(X)]=\psi(a/k)/k+\log(b) .
$$

（おしまい）