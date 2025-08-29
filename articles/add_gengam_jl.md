---
title: "Distributions.jl に自分で定義した分布を追加する；一般化ガンマ分布を例に"
emoji: "⛑️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [Julia]
published: false
---


次の密度関数 $f(x;a,b,k)$ を持つ確率分布を一般化ガンマ分布（[Generalized gamma distribution - wikipedia](https://en.wikipedia.org/wiki/Generalized_gamma_distribution)）という.

$$
f(x; a, b, k) = \frac{(k/b^a) x^{a-1} \exp(-(x/b)^k)}{\Gamma(a/k)}, \quad x > 0.
$$

$(a, b, k)$ はすべて正の値を取るパラメータである．Distributions.jl はいい感じのパッケージだが，幸か不幸か，一般化ガンマ分布は用意されていないようなので （ [Distributions.jl; Univariate Continuous Distributions](https://juliastats.org/Distributions.jl/stable/univariate/#Continuous-Distributions) ）これを実装してみた．具体的には，次のようにした．

あとは答えあわせのためにやったことをメモしておく．

パラメータ $(a, b, k)$ の値をいくつか選んでプロットしたものが下の図である.

確率変数 $X$ が一般化ガンマ分布に従うことを $X \sim \mathrm{GG}(a, b, k)$ と書くことにする.

## 付録：期待値と対数の期待値

$X \sim \mathrm{GG}(a, 1, k)$ の場合の期待値 $E[X]$ を定義通り計算すると次のようになる．

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

いま, $X$ の自然対数をとった確率変数 $\log X$ の期待値 $E[\log X]$ を求めたい（そういうのを求めたい状況がたまにある）.

その準備としてガンマ関数の定義を確認する．

$$
\Gamma(z) = \int^{\infty}_{0} t^{z-1}\exp(-t) \, dt
$$

である. $u = \log(t)$ と変数変換すると， $du = dt/t$ より, ガンマ関数は，次のように書ける.

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

$E[\log(bX)]=E[\log(X)]+\log(b)$ なので, $X \sim \mathrm{GG}(a, b, k)$ のときは $E[\log(X)]=\psi(a/k)/k+\log(b)$.

