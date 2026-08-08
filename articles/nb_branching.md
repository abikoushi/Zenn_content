---
title: "離散時間の分岐過程とその感染症疫学への応用"
emoji: "🪾"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 確率過程]
published: false
---

## このノートについて

このノートは [西浦 博; 小林 鉄郎; 安齋 麻美; 合原 一幸; ナタリー・リントン. 感染症流行を読み解く数理（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) 第5章の「MERSの伝播リスク」の節の行間を補うためのものです．

同書の (5.9) 式のシミュレーションによる確認と導出を行います． 同書 (5.9) 式がこのノートにおける (2) 式に対応します．

導出については，

Nishiura, H., Yan, P., Sleeman, C. K., & Mode, C. J. (2012). Estimating the transmission potential of supercritical processes based on the final size distribution of minor outbreaks. Journal of theoretical biology, 294, 48-55. https://www.sciencedirect.com/science/article/pii/S0022519311005625?via%3Dihub

も合わせて参考にしました．

また，確率母関数についての記述は [竹村彰通『現代数理統計学』（学術図書出版社）](https://www.gakujutsu.co.jp/product/978-4-7806-0860-1/) を参考にしました．

数値例は R によります．

## シミュレーションしてみる

1人の感染者が次の感染者を生み出し，その感染者がまた次の感染者を生み出し，……というような過程を考えたい．

最初の1人の感染者を第0世代の感染者と呼ぶことにしよう．第0世代の感染者から生み出された感染者を第1世代，第1世代の感染者から生み出された感染者を第2世代，……，第 $n-1$ 世代の感染者から生み出された感染者を第n世代の感染者と呼ぶことにする．

第 $n$ 世代の $i$ 番目の感染者が生み出す2次感染者数を $X_i(n)$ と置く．1人の感染者から生み出される感染者数 $X_i(n)$ を負の二項分布でモデル化する．1人の感染者から生み出される感染者数の平均が $R_0$ となるよう，ここでの負の二項分布の確率関数は次のようにパラメータ化される．

$$
\Pr(X_i(n)=x) = \frac{\Gamma(k+x)}{x!\Gamma(k)}\left(\frac{R_0}{R_0+k}\right)^x \left(1+\frac{R_0}{k}\right)^{-k} \tag{1}
$$

ここでの $\Gamma(x)$ はガンマ関数である．この分布の平均は $R_0$，分散は $R_0+R_0^2/k$ となる

R の `*nbinom` 関数を使う際は，次のように置くと (1) の確率関数と一致する．

$$
\texttt{size} \leftarrow k, \quad \texttt{prob} \leftarrow k/(k+R_0)
$$

イメージをつかむため，まずシミュレーションをやってみる．

結論を先取りする形で書くが，第 $\infty$ 世代までの総感染者数 $Y = \sum_{n=0}^{\infty}\sum_{i=1}X_i(n)$ の確率関数は次のようになる．

$$
\Pr(Y=y)=\begin{cases}
\frac{1}{(1+R_0/k)^k} & y=1\\
\frac{\prod_{j=0}^{y-2}(j/k+y)}{y!}\left(\frac{k}{R_0+k}\right)^{ky} \left(\frac{R_0k}{R_0+k}\right)^{y-1} & y\ge2
\end{cases} \tag{2}
$$

## 導出してみる

### 準備：確率母関数

確率関数 $f_X(x)$ を持つ確率変数 $X$ について，次のように定義される関数 $G_X(x)$ を確率母関数と呼ぶ．

$$
G_X(z) = E[z^X] = \sum_{x=0}^{\infty} z^x f_X(x) .
$$

ただし，ここでは $|z| \le 1$ とした．確率母関数 $G_X(z)$ が与えられれば，

$$
\frac{d^x}{dz^x}G_X(x) \mid_{z=0} = G^{(x)}_X(0) = (x!)f_X(x)
$$

より，

$$
f_X(x) = G^{(x)}_X(0)/x! 
$$

と，$X=x$ となる確率を求めることができる．

また，確率母関数をモーメントを求めるのに使うこともできる．たとえば，

$$
\begin{aligned}
G'(1) &= \sum_{x=0}^{\infty} x f_X(x) = E[X],\\
G''(1) &= \sum_{x=0}^{\infty} x(x-1) f_X(x) = E[X(X-1)].
\end{aligned}
$$

### 準備：負の二項分布

まず (1) の負の二項分布が確率関数の定義を満たすよう，総和が 1 となっていることを確かめよう．$p = k/(k+R_0)$, $q=1-p=R_0/(k+R_0)$ と置く． $p^{-k}=(1-q)^{-k}$ を $q$ の無限級数としてテイラー展開すると，

$$
\begin{aligned}
p^{-k}&=(1-q)^{-k}\\
&=1+kq+\frac{k(k+1)}{2!}q^2 + \frac{k(k+1)(k+2)}{3!} +\cdots \\
&=\sum_{x=0}^{\infty} \frac{\Gamma(k+x)}{\Gamma(k)x!} q^j.
\end{aligned}\tag{3}
$$

両辺に $p^{k}$ を掛けると (1) の総和が 1 となることがわかる．

次に負の二項分布の確率母関数を求める． (3) の関係より，次が成り立つ．

$$
\sum_{x=0}^{\infty} \frac{\Gamma(k+x)}{\Gamma(k)x!}(qz)^x = (1-qz)^{-k}
$$

これより，(1) の負の二項分布に従う確率変数 $X$ についての確率母関数は，

$$
\begin{aligned}
G_X(z)&=\sum_{x=0}^{\infty} \frac{\Gamma(k+x)}{\Gamma(k)x!} z^x q^x p^k\\
&=\left(\frac{p}{1-qz}\right)^k\\
&=\left(\frac{1-qz-q+q}{p}\right)^{-k}\\
&=\left(1-(z-1)\frac{q}{p}\right)^{-k}\\
&=\left(1-(z-1)\frac{R_0}{k}\right)^{-k}.
\end{aligned}
$$


### 総感染者数の分布


第 $n$ 世代の感染者数の合計を $X(n)$ と書くことにする．第 $n$ 世代までの感染者数の和を $Y_n = \sum_{j=1}^n\sum_{i=1}^{X(n)}X_i(j)$ と書くことにする．
(1) の確率関数を $f_X(x)$ と置く．



第0世代から第1世代までの総感染者数の合計は $Y_1 = 1+X$ であり，その確率母関数は，

$$
G_{Y_1}(z) = E[z^{1+X}f_X(x)] = z G_X(z).
$$



$$
G_{Y_2}(z) = E[z^{Y_1+\sum_{i=1}^{X(2)}X_i(2)}f_X(x)] = z G_X(z).
$$