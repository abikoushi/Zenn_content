---
title: "スピアマンとケンドールの順位相関係数"
emoji: "📑"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: true
---

## このノートについて

順序統計量を使う相関係数として有名な2つ，**スピアマン（Spearman）の順位相関係数** と **ケンドール（Kendall）の順位相関係数** を紹介する．

といってもほとんど計算メモのようなものである．

順位相関係数の標準的な定義としては，主に [竹内啓『数理統計学―データ解析の方法』（東洋経済新報社）](https://books.google.co.jp/books/about/数理統計学.html?id=fSo8DwAAQBAJ&redir_esc=y) を参照した．

タイ（同順位）があるデータのときはちょっと難しくなるのでここでは扱わない．

また，このノートで述べているものは標本相関係数だけであり，確率論までは扱わない．


## ピアソンの積率相関係数からスピアマンの順位相関係数へ


単に相関係数と言うときは **ピアソン（Pearson）の積率相関係数** を指すことが多い．組のデータ $(x_i，y_i)$，($i=1,\ldots n$) が得られたとき，ピアソンの積率相関係数 $r$ は次の式で求められる．

$$
r = \frac{\sum_{i=1}^n (x_i - \bar x)(y_i - \bar y)}{\sqrt{\sum_{i=1}^n (x_i - \bar x)^2}\sqrt{\sum_{i=1}^n (y_i - \bar y)^2}} \tag{P}
$$

ただしここでは $\bar x= \frac{1}{n}\sum_{i=1}^n x_i$ とした.

一方で，組のデータの順位 $(x_i，y_i)$，($i=1,\ldots n$) が得られたとき，スピアマンの順位相関係数 $\rho$ は次の式で求められる．

$$
\rho = 1-\frac{6\sum_{i=1}^n (x_i - y_i)^2}{n(n^2-1)} \tag{S}
$$

ピアソンの積率相関係数が色々な相関係数の基本のように思えるので，まずピアソンの積率相関係数の式の意味を一度考えてみる．

文字が多いと大変なので，仮に $\bar x$ と $\bar y$ がどちらも 0 としてみる．

$$
r = \frac{\sum_{i=1}^n x_iy_i}{\sqrt{\sum_{i=1}^n x_i^2}\sqrt{\sum_{i=1}^n y_i^2}}. \tag{P0}
$$

これは2つのベクトル $x=(x_1,\ldots x_n)'$, $y=(y_1,\ldots y_n)'$ の**内積** 

$$
\langle x,y\rangle = \sum_{i=1}^n x_iy_i
$$ 

をそれぞれの**長さ**

$$
\|x\|=\sqrt{\sum_{i=1}^n x_i^2}, \quad \|y\|=\sqrt{\sum_{i=1}^n y_i^2}
$$

で割ったもの，つまり2つのベクトルのなす角のコサインになっている．そのため (P0) の $r$ を **コサイン類似度（cosine similalrity）** と呼ぶこともある．

$\bar x = 0$，$\bar y=0$ は $x_i$，$y_i$ からそれぞれの標本平均を引いたものを改めて $x_i$，$y_i$ とみれば達成されるので，相関係数は平均からのズレ具合についてのコサイン類似度である．

コサインであるから相関係数は $-1\le r \le 1$ である．

ここでは確率論についてはほとんど触れないので，線形代数の教科書と数理統計の教科書を見比べながら勉強する場合のために，模式的に対応を書いておく．

|数ベクトルの場合|確率変数の場合|
|--|--|
|長さ|標準偏差|
|内積|共分散|
|なす角 $\theta$ の $\cos(\theta)$|相関係数|

例えば，内積についてのコーシー＝シュワルツの不等式に関しては，

$$
\begin{aligned}
\langle x,y\rangle ^{2} &\leq \langle x,x\rangle \cdot \langle y,y\rangle \\
& =\|x\|^2 \cdot \|y\|^2
\end{aligned}
$$

$x$, $y$ を確率変数 $X$, $Y$ として，内積を共分散 $\mathrm{Cov}(X,Y)$ と分散 $\mathrm{Var}(X)=\mathrm{Cov}(X,X)$ に置き換え，

$$
\begin{aligned}
\mathrm{Cov}(X,Y)^2 &\leq \mathrm{Cov}(X,X) \cdot \mathrm{Cov}(Y,Y)\\
&=\mathrm{Var}(X) \cdot \mathrm{Var}(Y)
\end{aligned}
$$

としたものも成り立ち，これもコーシー＝シュワルツの不等式と呼ぶ．

さて， (S)の式のスピアマンの順位相関係数には2つのベクトルの距離の2乗

$$
d^2(x,y) = \sum_{i=1}^n (x_i - y_i)^2
$$

が現れている．

２つのベクトル $x$, $y$ の距離 $d(x,y)$ は次のようにも書ける．

$$
\begin{aligned}
d(x,y) &= \sqrt{\sum_{i=1}^n (x_i-y_i)^2}\\
&=\sqrt{\left(\sum_{i=1}^n x_i^2\right)+\left(\sum_{i=1}^ny_i^2\right) - 2\left(\sum_{i=1}^nx_iy_i\right)}\\
&=\sqrt{\| x \|^2 + \|y\|^2 - 2\langle x, y\rangle}
\end{aligned}
$$

ルートの中の第1項，第2項は0以上であり，内積 $\langle x,y\rangle$ を含む項の符号がマイナスになっているので，内積が大きいほど距離は小さく，内積が小さいほど距離は大きくなることがわかる．

内積を，次のように距離が明示的に現れる形で書くこともできる．

$$
\langle x, y\rangle =\frac{1}{2} \left(  \| x \|^2 + \|y\|^2  -d^2(x,y) \right)
$$

これを使って，相関係数 $r$ を距離が明示的に現れる形で書いてみる．こんなふうだ．

$$
\begin{aligned}
r &= \frac{\left(\sum_{i=1}^n x_iy_i \right) - n\bar{x}\bar{y}}{(\left\{ \sum_{i=1}^n x_i^2 \right\} - n\bar{x}^2)^{1/2} \cdot (\left\{ \sum_{i=1}^n y_i^2\right\} - n\bar{y}^2)^{1/2}}\\
&=\frac{\langle x, y\rangle - n\bar{x}\bar{y}}{(\|x \|^2- n\bar{x}^2)^{1/2} \cdot (\| y\|^2 - n\bar{y}^2)^{1/2}}\\
&=\frac{\frac{1}{2} \left\{\|x\|^2 + \|y\|^2 - d^2(x,y) \right\}- n\bar{x}\bar{y}}{(\|x \|^2- n\bar{x}^2)^{1/2} \cdot (\| y\|^2 - n\bar{y}^2)^{1/2}}
\end{aligned}
$$

いま，$x_i$，$y_i$ はどちらも1から $n$ までの重複のない自然数とし，データに値の小さい順に振られた順位を表すとする．このとき，次の2つが成り立つ．

$$
\begin{aligned}
n\bar{x}= \sum_{i=1}^n x_i & =\sum_{i=1}^n y_i = n\bar{y} ,\\
\|x\|^2 = \sum_{i=1}^n x_i^2 &=\sum_{i=1}^n y_i^2 = \|y\|^2.
\end{aligned}
$$

そのため，

$$
\begin{aligned}
r &= \frac{\frac{1}{2} \left\{  2\| x \|^2 - d^2(x,y) \right\}- n\bar{x}^2}{\|x \|^2- n\bar{x}^2}\\
&= 1 - \frac{\frac{1}{2}d^2(x,y)}{\|x \|^2- n\bar{x}^2}
\end{aligned}
$$

となる．さらに右辺第2項の分母について整理する．

$$
\begin{aligned}
\text{(分母)} &= \left(\sum_{i=1}^n x_i^2\right) - \frac{1}{n}\left(\sum_{i=1}x_i\right)^2\\
&= \left(\sum_{i=1}^n i^2\right) - \frac{1}{n}\left(\sum_{i=1}i\right)^2\\
&= \frac{1}{6} n (n+1)(2n+1) - \frac{1}{n}\left(\frac{1}{2}n(n+1)\right)^2\\
&=\frac{1}{12} ( n(n+1) [ 2(2n+1 )- 3(n+1)] )\\
&=\frac{1}{12} n(n^2-1)
\end{aligned}
$$

結果，$(x_i，y_i)$ が順位のときは，

$$
r = \rho = 1-\frac{6\sum_{i=1}^n (x_i - y_i)^2}{n(n^2-1)}
$$

となり (S) と一致することがわかった．すなわち，スピアマンの順序相関係数は，単に $(x_i，y_i)$ を順位に変換してからピアソンの積率相関係数を求めたものだとわかった．

一般に，$0 \le d^2(x，y) \le D_{\max}$ ならば

$$
-1 \le 1 - \frac{d^2(x，y)}{D_{\max}/2} \le 1 \tag{I}
$$

である． 符号を反転させた $\frac{d^2(x，y)}{D_{\max}/2}-1$ も $[-1,1]$ の範囲に収まるが，距離が大きいほど相関係数は小さくなってほしいので相関係数としては $1 - \frac{d^2(x，y)}{D_{\max}/2}$ のほうを採用する．

$d^2(x，y)$ の取りうる最大値 $D_{\max}$ は

$$
D_{\max} = \sum_{i=1}^n(i-(n-i+1))^2 = n(n^2-1)/3
$$

であるから，ここからも $\rho$ が (I)の不等式を満たすことがわかる．

## スピアマンの順位相関係数からケンドールの順位相関係数へ

順序統計量を使う相関係数として有名なものにもう一つ，ケンドールの順位相関係数がある．

スピアマンの順位相関係数では距離の2乗として差の2乗の総和を採用していた．

統計モデルを理解する上では，データが同じでも距離の測り方をいろいろ変えてみることができるという発想が大事になるように思う．

そこで，$x$ を $y$ と一致させるために必要な入れ替えの回数を距離の2乗としてみよう．$K$ 回の入れ替えが必要なとき距離 $d^2(x,y)=K$ である．

たとえば，$x=(1,3,2)'$ , $y=(1,2,3)'$ なら，2つめと3つめを1回入れ替えると一致するので $d^2(x,y)=1$ である.

一般に，もれなくダブりなく何回の入れ替えが必要か数えるのはどうしたらいだろうか．

次の表を例に考える．

|x|y|
|--|--|
|3|4|
|1|2|
|4|3|
|2|1|

まず $y$ をキーにソートする．

|x|y|
|--|--|
|2|1|
|1|2|
|4|3|
|3|4|


表の $i$ 行目より下にあるもののうち，$x_i$ より大きいものは入れ替えが必要がない．この数 $c_i$ の合計を $S=\sum_{i=1}^n c_i$ とする．

先の例では，$S$ は次の表の太字の部分の数で，$c_1=2$, $c_2=2$, $c_3=0$, $c_4=0$ である．

|1行目に注目|2行目に注目|3行目に注目|4行目に注目|
|--|--|--|--|
|2||||
|1|1|||
|**4**|**4**|4||
|**3**|**3**|3|3|

最大の入れ替え回数は二項係数を使って $\binom{n}{2}=n(n-1)/2$ と書けるので，

$$
K=n(n-1)/2-S
$$

である．先の例では $K=6-4=2$ である．

$0 \le K \le n(n-1)/2$ より(I)の関係を使って，

$$
\tau = 1 - \frac{K}{n(n-1)/4}
$$

とすると $-1 \le \tau \le 1$ を満たす． $\tau$ をケンドールの順位相関係数と呼ぶ．

ケンドールの順位相関係数は，スピアマンの順位相関係数において，距離の2乗を「2つの順位を一致させるために必要な入れ替えの回数」（この距離にはたぶんなにか名前がついていると思うが見つからなかった……と思ったけど追記．『ケンドール距離』でいいらしいです：[神嶌敏弘 - 順序の距離と確率モデル](https://www.jstage.jst.go.jp/article/jsaisigtwo/2009/DMSM-A902/2009_07/_article/-char/ja/) ）に変えたものである．

## R による計算例

R の `cor` 関数は `method = "spearman"` や `method = "kendall"` も実装されている（デフォルトは `method = "pearson"`）．

あまりおもしろい例が作れなかったが，こんな感じで答え合わせができる．

```r
set.seed(0705)
n=11
x = rnorm(n)
y = rnorm(n)

rho0 = cor(x, y, method = "spearman")
D2 = sum((rank(x)-rank(y))^2)
D2max = n*(n+1)*(n-1)/3
rho = 1-D2/(D2max/2)
cat("our implementation", rho, ", stats::cor", rho0)

###

tau0 = cor(x, y, method = "kendall")
R = rank(x[order(y)])
S = sum(sapply(1:(n-1),function(i)sum(R[(i+1):n] > R[i])))
K = choose(n,2)-S
Kmax = choose(n,2)
tau = 1-K/(Kmax/2)

cat("our implementation", tau, ", stats::cor", tau0)
```

```r
cat("our implementation", rho, ", stats::cor", rho0)
our implementation -0.06363636 , stats::cor -0.06363636

cat("our implementation", tau, ", stats::cor", tau0)
our implementation -0.01818182 , stats::cor -0.01818182
```

おしまい．