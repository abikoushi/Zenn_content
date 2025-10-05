---
title: "二項分布とベータ分布の関係+α"
emoji: "🗼"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: false
---

## まえおき

[ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) のいとこみたいな記事です．ただし [ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) を先に読んでおく必要はありません．

結論を先取りする形で確率分布どうしの関係を表にまとめると，こんな感じになります．

|  | 上限が自由 | 上限が固定 |
|---|---|---|
| イベントの回数 | ポアソン分布 | 二項分布 |
| イベントの待ち時間 | ガンマ分布| ベータ分布 |

しかしこの表を見て意味がわかる人はすでにわかってる人だと思うので，以下の本文に続きます．


## 動機のための例：信頼区間の紹介

二項分布のモデル $X \sim \mathrm{Binom}(n,p)$ を考える．$p$ の信頼区間を作るとき，次の関係を用いることができる．

$$
\sum_{k=x}^n \binom{n}{k} p^k (1-p)^{n-k}
= \frac{1}{\Beta(x+1)} \int_0^{p} u^{x-1} (1-u)^{n-x-1}\, du.
\tag{1}
$$



## 数値計算してみる

(1) は次のように捉えるとわかりやすい．区間 $[0,1]$ の一様乱数を $n$ 個生成したとき，

$$
\begin{aligned}
P&(\text{ {\it t} 以下のイベントが {\it k} 回以上発生する})\\
&＝P(\text{小さいほうから {\it k} 番目のイベントが {\it t} 以下}).
\end{aligned}
$$



## 手計算してみる

(1) の左辺を $G(p)$ とおく．

$$
G(p) = \sum_{k=m}^n \binom{n}{k} p^k (1-p)^{n-k}
$$

その導関数は

$$
\begin{aligned}
G'(p) &= \sum_{k=m}^n \binom{n}{k} 
   \left[ k p^{k-1}(1-p)^{n-k} - p^k (n-k)(1-p)^{n-k-1} \right] \\
&= \sum_{k=m}^n \frac{n!}{(n-k)!(k-1)!} p^{k-1}(1-p)^{n-k} - \sum_{k=m}^n \frac{n!}{(n-k-1)!k!}  p^k (1-p)^{n-k-1}
\end{aligned}
$$

ここで第2項の添字を $j = k+1$ と置き換えると，

$$
G'(p)= \sum_{k=m}^n \frac{n!}{(n-k)!(k-1)!} 
   p^{k-1}(1-p)^{n-k}
 - \sum_{j=m+1}^{n+1} \frac{n!}{(n-j)!(j-1)!} 
   p^{j-1}(1-p)^{n-j}
$$


結果，$k=m$ の項だけが残り，

$$
G'(p) = \frac{n!}{(n-m)!(m-1)!} p^{m-1}(1-p)^{n-m}.
$$

したがって，

$$
G'(p) = \frac{1}{B(m, n-m+1)} p^{m-1}(1-p)^{n-m},
$$

となる．これを $p$ で積分すれば，

$$
\begin{aligned}
G_m(p) &= \int_0^p G'(u)\, du\\
&= \frac{1}{B(m, n-m+1)} \int_0^p 
     u^{m-1}(1-u)^{n-m}\, du,
\end{aligned}
$$

であり，これはベータ分布の分布関数である．(1) がわかった．

[竹村彰通『新装改訂版 現代数理統計学』学術図書出版社](https://www.gakujutsu.co.jp/product/978-4-7806-0860-1/) の「4.7 」

## プラスアルファの関係

実は，ここまで述べたことは [ポアソン分布とガンマ分布の関係](https://zenn.dev/abe2/articles/rel_poisson_and_gamma) とほぼ同じ構成だった．それを思うと，二項分布とポアソン分布，ガンマ分布とベータ分布が互いにそれぞれ関係が深いように見える．そこで二項分布とポアソン分布，ガンマ分布とベータ分布の関係についても触れておく．


### ポアソン分布と二項分布の関係

$X_i$ ($i=1,2$) がそれぞれ独立にパラメータ $\lambda_i$ のポアソン分布 $\mathrm{Pois}(\lambda_i)$ に従うとする．

$S=X_1+X_2$ とおき，和 $S$ が適当な値 $n$ に固定されたときの $X_1$ の分布を求める．

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

$S=X_1+X_2$ とおき，和が $S=1$  に固定されたときの $X_1$ の分布を求める．

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

この節はおまけみたいなもの．統計学の教科書には上述のベータ分布ではなくF分布を使って母比率の信頼区間を表示している場合がある（例えば [野田一雄・宮岡悦良『入門・演習 数理統計』共立出版](https://www.kyoritsu-pub.co.jp/book/b10011485.html) ）．そこでベータ分布とF分布の関係についても述べる．

$X \sim \mathrm{Beta}(a,b)$ のとき， $Y = X/(1-X)$ とすると, $Y$ の分布関数は，$X = Y/(1+Y)$, $1-X = 1/(1+Y)$, $dx = (1+Y)^{-2} dy$ より

$$
\begin{aligned}
&P(Y \leq u) 
   = \int_0^u \frac{1}{B(a,b)} \left(\frac{y}{1+y}\right)^{a-1}
       \left(\frac{1}{1+y}\right)^{b-1} \frac{dy}{(1+y)^2} \\
&= \int_0^u \frac{1}{B(a,b)} y^{a-1}(1+y)^{-(a+b)} \, dy 
\end{aligned}
$$

これはベータプライム分布（逆ベータ分布，第2種ベータ分布とも呼ばれる；[Wikipedia - Beta prime distribution](https://en.wikipedia.org/wiki/Beta_prime_distribution)）の分布関数である．

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

これはF分布（[Wikipedia - F-distribution](https://en.wikipedia.org/wiki/F-distribution)）の分布関数である．

これは私見だが，F分布を使って書けるようにしているのは巻末のF分布表を引けるように工夫した結果と推測できる．コンピューターを使って計算することが主な現代においては，あまり重要視しなくていいように思う．
