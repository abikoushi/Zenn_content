---
title: "二項分布とベータ分布の関係+α"
emoji: "🗼"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: false
---

## まえおき

[ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) のいとこみたいな記事です．ただし [ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) を先に読んでおく必要はありません．

結論を先取りする形で確率分布どうしの関係を表にまとめるとこんな感じになります．

|  | 上限が自由 | 上限が固定 |
|---|---|---|
| イベントの回数 | ポアソン分布 | 二項分布 |
| イベントの待ち時間 | ガンマ分布| ベータ分布 |

しかしこの表を見て意味がわかる人はすでにわかってる人だと思うので，以下の本文に続きます．

この文書を書くのに使用した R のコードは末尾にまとめてあります．

## 動機のための例：信頼区間の紹介

二項分布のモデル $X \sim \mathrm{Binom}(n,p)$ を考える．パラメータ $p$ についての信頼区間を作るとき，次の関係を用いることができる．

$$
\sum_{k=x}^n \binom{n}{k} p^k (1-p)^{n-k}
= \frac{1}{\Beta(x,n-k+1)} \int_0^{p} u^{x-1} (1-u)^{n-x}\, du.
\tag{1}
$$

左辺を $x$ の関数とみると二項分布の補分布関数であり，右辺を $p$ の関数とみるとベータ分布の分布関数である．

次のように適当な数を選んで R で確認してみよう．

```r
n = 10L
k = 4L
p = 0.6
```

```r
> pbinom(k-1, n, p, lower.tail = FALSE)
[1] 0.9452381
> pbeta(p, shape1 = k, shape2 = n-k+1)
[1] 0.9452381
```

(1)より次の(2)も成り立つ．

$$
\sum_{k=1}^x \binom{n}{k} p^k (1-p)^{n-k}
= \frac{1}{\Beta(x+1,n-x)} \int_{p}^{1} u^{x} (1-u)^{n-x-1}\, du.
\tag{2}
$$

$Q(\alpha;a,b)$ をベータ分布 $\mathrm{Beta}(a,b)$ の分布関数の逆関数とすると，(1)と(2)を用いて信頼区間 $C(X)$ を次のように作ることができる．

$$
C(X) = [Q(\alpha/2; X+1,n-X), \ Q(1-\alpha/2;X,n-X+1)] .
$$

この信頼区間は R では `binom.test` という関数で実装されている．

```r
> res = binom.test(k, n, p=p)
> res$conf.int
[1] 0.1215523 0.7376219
attr(,"conf.level")
[1] 0.95
> c(qbeta(0.025, shape1 = k, shape2 = n-k+1),
+   qbeta(0.975, shape1 = k+1, shape2 = n-k))
[1] 0.1215523 0.7376219
```

加えて，ベイズ信頼区間も考えてみる．

$p$ の事前分布としてベータ分布 $\mathrm{Beta}(a,b)$ を用いると事後分布はベータ分布 $\mathrm{Beta}(x+a,n-x+b)$ であり，ベイズ信頼区間 $B(X)$ は次のように得られる．

$$
B(X) = [Q(\alpha/2;X+a,n-x+b), \ Q(1-\alpha/2;X+a,n-X+b)] .
$$

$C(X)$ と $B(X)$ は互いに似ている．

信頼区間をプロットしてみるとベイズ版のほうが狭めになる．事前分布のパラメータは一様分布に相当する $a=1$, $b=1$ に設定した．

![](/images/rel_/pvfun_binom1.png)
*縦軸を $\alpha$ ， 横軸を信頼区間とした関数（信頼区間関数 a.k.a. P値関数）． k=4, n=10のとき*

標本サイズ $n$ を増やすとふたつの違いは小さくなる．

![](/images/rel_/pvfun_binom2.png)
*縦軸を $\alpha$ ， 横軸を信頼区間とした関数（信頼区間関数 a.k.a. P値関数）． k=40, n=100のとき*

このように出発点の考え方が違っても似た結果が得られるというのは，数理的なことを学ぶ上でおもしろいポイントの一つだと思う．

以上により，（1）の関係を理解したい動機が得られたことにして，次節のシミュレーションに進む．


## 数値計算してみる

(1) は次のように捉えるとわかりやすい．区間 $[0,1]$ の一様乱数を $n$ 個生成したとき，

$$
\begin{aligned}
P&(\text{ {\it p} 以下の乱数が {\it k} 回以上発生する})\\
&＝P(\text{小さいほうから {\it k} 番目の乱数の値が {\it p} 以下}).
\end{aligned}
$$

「$p$ 以下の乱数が $k$ 回以上発生する」と「小さいほうから $k$ 番目の乱数が $p$ 以下」は同値である．

そして、次のようにシミュレーションを書くと，「$p$ 以下の乱数の発生回数」と「小さいほうから $k$ 番目の乱数の値」がそれぞれ二項分布とベータ分布に従うことが確かめられる．

```r
n = 10L
k = 4L
p = 0.6
simorderstat_unif <- function(iter, k, p, n){
  res_p = numeric(iter)
  res_k = integer(iter)
  for(it in seq_len(iter)){
    u = dqrng::dqrunif(n) #faster than stats::runif
    res_p[it] = sort(u)[k]
    res_k[it] = sum(u<p)
  }
  return(data.frame(waitingtime=res_p, counts=res_k))
}
```

![](/images/rel_binom_beta/simcount_binom.png)
*青いマーカーは二項分布の確率関数*

![](/images/rel_binom_beta/simtime_beta.png)
*青い点線はベータ分布の密度関数*

ここでは，$k$ 番目に小さい値を $k$ 番目までの待ち時間（waiting time）ということにした．


## 手計算してみる

$p$ 以下のイベントの発生回数が二項分布になるのはいいとして，小さいほうから $k$ 番目の値がベータ分布になるのは意外な感じがするかもしれない．そこで(1)の証明もみておこう．

(1) の左辺を $G(p)$ とおく．

$$
G(p) = \sum_{k=m}^n \binom{n}{k} p^k (1-p)^{n-k}
$$

その導関数は

$$
\begin{aligned}
G'(p) &= \sum_{k=x}^n \binom{n}{k} 
   \left[ k p^{k-1}(1-p)^{n-k} - p^k (n-k)(1-p)^{n-k-1} \right] \\
&= \left[\sum_{k=x}^n \frac{n!}{(n-k)!(k-1)!} p^{k-1}(1-p)^{n-k}\right] - \left[\sum_{k=x}^n \frac{n!}{(n-k-1)!k!}  p^k (1-p)^{n-k-1} \right]
\end{aligned}
$$

ここで第2項の添字を $j = k+1$ と置き換えると，

$$
G'(p)= \left[ \sum_{k=x}^n \frac{n!}{(n-k)!(k-1)!} p^{k-1}(1-p)^{n-k}\right]
 - \left[ \sum_{j=x+1}^{n+1} \frac{n!}{(n-j)!(j-1)!} p^{j-1}(1-p)^{n-j}\right]
$$


結果，$k=x$ の項だけが残り，

$$
G'(p) = \frac{n!}{(n-x)!(x-1)!} p^{x-1}(1-p)^{n-x}.
$$

したがって，

$$
G'(p) = \frac{1}{B(m, n-m+1)} p^{m-1}(1-p)^{n-m},
$$

となる．これを積分すれば，

$$
\begin{aligned}
G(p) &= \int_0^p G'(u)\, du\\
&= \frac{1}{B(m, n-m+1)} \int_0^p 
     u^{m-1}(1-u)^{n-m}\, du,
\end{aligned}
$$

であり，これはベータ分布の分布関数である．(1) がわかった．

これと同様のことが [竹村彰通『新装改訂版 現代数理統計学』学術図書出版社](https://www.gakujutsu.co.jp/product/978-4-7806-0860-1/) の「4.7 順序統計量と経験分布関数」にもう少し一般化した形で書かれているので，この節を書くにあたって参考にした．


## プラスアルファの関係

実は，ここまで述べたことは [ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) とほぼ同じ構成であった．それを思うと，二項分布とポアソン分布，ガンマ分布とベータ分布が互いにそれぞれ関係が深いように見える．そこで二項分布とポアソン分布，ガンマ分布とベータ分布の関係についても触れておく．


### ポアソン分布と二項分布の関係

$X_i$ ($i=1,2$) がそれぞれ独立にパラメータ $\lambda_i$ のポアソン分布 $\mathrm{Pois}(\lambda_i)$ に従うとする．

$S=X_1+X_2$ とおき，再生性 $S \sim \mathrm{Pois}(\lambda_1+\lambda_2)$ を使って，和 $S$ が適当な値 $n$ に固定されたときの $X_1$ の分布を求める．

$$
\begin{aligned}
&P(X_1=x | S = n) = \frac{P(X_1=x) P(X_2 = n-x)}{P(S=n)}\\
&= \left. \left[ \frac{\lambda_1^x e^{-\lambda_1}}{x!} \cdot \frac{\lambda_2^{n-x} e^{-\lambda_2}}{(n-x)!} \right]  \middle/ \left[\frac{\{\lambda_1+\lambda_2\}^n e^{-(\lambda_1+\lambda_2)}}{n!} \right]\right.\\
&= \frac{n!}{x! (n-x)!} 
   \left( \frac{\lambda_1}{ \lambda_1+\lambda_2} \right)^{x_j}
   \left( \frac{\lambda_2}{\lambda_1+\lambda_2} \right)^{n-x_j}
\end{aligned}
$$

これは二項分布の確率関数である．


### ガンマ分布とベータ分布の関係

$X_i$ ($i=1,2$) がそれぞれ独立にパラメータ ($a_i$, $b$) のガンマ分布 $\mathrm{Gamma}(a_i, b)$ に従うとする．

$S=X_1+X_2$ とおき，再生性 $S \sim \mathrm{Gamma}(a_1+a_2, \, b)$ を使って，和が $S=1$  と固定されたときの $X_1$ の分布を求める．

$$
\begin{aligned}
&P\!\left(X_1 = x_1 | S = 1 \right)
=  \frac{ P\!\left(X_1 = x_1 \right) \; P\!\left( X_2 = 1 - x_1 \right) }{ P( S = 1 ) }\\
&= \left. \left[ \dfrac{b^{a_1}}{\Gamma(a_1)} x_1^{a_1 - 1} e^{-b x_1}
        \cdot \dfrac{b^{a_1+a_2}}{\Gamma(a_1+a_2)}\right] \middle/ \left[
         (1-x_1)^{a_1+a_2 - 1} e^{ -b (1-x_1)} 
       { \dfrac{b^{a_1+a_2}}{\Gamma\!\left(a_1+a_2\right)} }\right]\right.\\
&= \frac{\Gamma(a_1 + a_2)}{\Gamma(a_1)\,\Gamma(a_2)}
   x_1^{a_1 - 1} (1-x_1)^{ a_2 - 1}
\end{aligned}
$$

これはベータ分布の密度関数である．

これまでに述べた関係を表にまとめると次のようになる．

|  | 上限が自由 | 上限が固定 |
|---|---|---|
| イベントの回数 | ポアソン分布 | 二項分布 |
| イベントの待ち時間 | ガンマ分布| ベータ分布 |

これは冒頭の表の再掲である．


### ベータ分布とF分布の関係

この節はおまけみたいなもの．統計学の教科書には上述のベータ分布ではなくF分布を使って母比率の信頼区間を表示している場合がある（例えば [野田一雄・宮岡悦良『入門・演習 数理統計』共立出版](https://www.kyoritsu-pub.co.jp/book/b10011485.html) の6章の章末問題）．そこでベータ分布とF分布の関係についても少し述べる．

これは私見だが，F分布を使っているのは巻末のF分布表を引いて計算できるようにしたためと思われる．コンピュータを使って計算することが主な現代においてはあまり重要視しなくてもいいだろう．

$X \sim \mathrm{Beta}(a,b)$ のとき， $Y = X/(1-X)$ とすると, $Y$ の分布関数は，$X = Y/(1+Y)$, $1-X = 1/(1+Y)$, $dx = (1+Y)^{-2} dy$ より

$$
\begin{aligned}
&P(Y \leq u) 
   = \int_0^u \frac{1}{B(a,b)} \left(\frac{y}{1+y}\right)^{a-1}
       \left(\frac{1}{1+y}\right)^{b-1} \frac{dy}{(1+y)^2} \\
&= \int_0^u \frac{1}{B(a,b)} y^{a-1}(1+y)^{-(a+b)} \, dy 
\end{aligned}
$$

これはベータプライム分布（[Wikipedia: Beta prime distribution](https://en.wikipedia.org/wiki/Beta_prime_distribution)）の分布関数である．逆ベータ分布，第2種ベータ分布とも呼ばれる．

さらに $Z = (a/b) \, Y$ とおくと, $Y = (b/a) \, Z$, $dy = (b/a) \, dz$ より，$Z$ の分布関数は，

$$
\begin{aligned}
P(Z \leq u) &= \int_0^u \frac{1}{B(a,b)} \left(\tfrac{b}{a}z\right)^{a-1}
       \left(1+\frac{b}{a}z\right)^{-(a+b)} \frac{b}{a} \, dz \\
&= \int_0^u \frac{1}{B(a,b)} \left(\frac{b}{a}\right)^a 
       z^{a-1}\left(1+\frac{b}{a}z\right)^{-(a+b)} dz.
\end{aligned}
$$

パラメータを $a = d_1/2$,  $b = d_2/2$ とおくと, 

$$
P(Z \leq u)   = \int_0^u \frac{1}{B(d_1/2,d_2/2)} 
     \left(\frac{d_1}{d_2}\right)^{d_1/2}
     z^{d_1/2-1}
     \left(1+\frac{d_1}{d_2}z\right)^{-(d_1+d_2)/2} dz.
$$

これはF分布（[Wikipedia: F-distribution](https://en.wikipedia.org/wiki/F-distribution)）の分布関数である．


## R のコード

https://github.com/abikoushi/Zenn_content/blob/main/R/rel_binom_and_beta.R