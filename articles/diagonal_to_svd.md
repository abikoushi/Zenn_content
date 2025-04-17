---
title: "行列の対角化から特異値分解までの行間を埋めるためのノート"
emoji: "↘️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 数値計算]
published: true
---


行列の対角化の話は線形代数の教科書として指定されるような本にはだいたい載っていると思う．

一方で統計・データ分析などの分野でよく使われる [特異値分解](https://ja.wikipedia.org/wiki/特異値分解) という言葉は線形代数の教科書にはあまり出ていないと思う．

このノートでは少しだけそのギャップを埋める手助けを目指す．

## 手計算してみる

まず線形代数の教科書によく出てくる対角化の方を紹介するが，その前に準備として用語を確認しておく．

正方行列 $A$ に対して

$$
A\boldsymbol{x} = \alpha \boldsymbol{x}
$$

となるベクトル $\boldsymbol{x}$ を **固有ベクトル** (eigen vector)，$\alpha$ を **固有値** (eigen value) という．

上で対角化と書いていたのは次のようなことだ．

:::message
**定理 1** ：実対称行列 $A$ に対して，ある直交行列 $P$ があり，

$$
P'AP = P^{-1}AP =\left(\begin{matrix} \alpha_1 & 0 & 0 & 0\\
0 & \alpha_2 & 0 & 0\\
\vdots & \ddots & \ddots & \vdots\\
0 & 0 & 0 &  \alpha_r  \\
 \end{matrix}\right)
$$

とできる．ここで $\alpha_i$ は $A$ の固有値とした．
:::

定理 1 は線形代数の教科書によく書かれている．もし手元に線形代数の本があれば「固有値」や「対角化」や「標準形」の章を参照してみてほしい．具体的には少なくとも下記にリストした文献には記載があることを確認した．

- [水田義弘『理工系　線形代数』（サイエンス社）](https://www.saiensu.co.jp/search/?isbn=978-4-7819-0859-5&y=1997)
- [基礎数学研究会　編『線形代数』（東海大学出版会）](https://www.hanmoto.com/bd/isbn/9784486012627)
- [佐竹一郎『線形代数学（新装版）』（裳華房）](https://www.shokabo.co.jp/mybooks/ISBN978-4-7853-1316-6.htm)
- [長谷川浩司『線形代数［改訂版］』（日本評論社）](https://www.nippyo.co.jp/shop/book/6704.html)


定理 1 から次の定理 2 が言える．

:::message
**定理 2**：$m$ 行 $n$ 列の任意の実行列 $A$ に対して $m$ 次の直交行列 $P_1$ と $n$ 次の直交行列 $P_2$ があり，

$$
P_1 A P_2 = \left(\begin{matrix} \gamma_1 & 0 & 0 & 0\\
0 & \gamma_2 & 0 & 0\\
\vdots & \ddots & \ddots & \vdots\\
0 & 0 & 0 &  \gamma_r  \\
 \end{matrix}\right)
$$

とできる．ここで $\alpha_i = \gamma_i^2$ は $(A' A)$ の固有値とした．$\gamma_1, \ldots, \gamma_r$ を特異値と呼ぶ．$r$ は $(A' A)$ のランクである．
:::

定理 2 が言えると

$$
A  =P_1' \left(\begin{matrix} \gamma_1 & 0 & 0 & 0\\
0 & \gamma_2 & 0 & 0\\
\vdots & \ddots & \ddots & \vdots\\
0 & 0 & 0 &  \gamma_r  \\
 \end{matrix}\right) P_2'
$$

も言える．これは行列 $A$ の特異値分解である．


### 定理 1 から 定理 2 を導く

$A' A$ は実対称行列なので定理 1 から

$$
P_2' (A' A) P_2 = \left(\begin{matrix} \alpha_1 & 0 & 0 & 0\\
0 & \alpha_2 & 0 & 0\\
\vdots & \ddots & \ddots & \vdots\\
0 & 0 & 0 &  \alpha_r  \\
\end{matrix}\right)
$$

とできる．

方針としては，$P_2' (A' A) P_2 =  (A P_2)'  A P_2$ なので $(A P_2)$ を 1 つの行列とみればいい．

ただし直交行列になるように，$\gamma_i$ で割っておく．

つまり，

$$
(A P_2) = (\boldsymbol{t}_1, \ldots, \boldsymbol{t}_r)
$$

としたとき，$P_1$ を

$$
P_1 = (\boldsymbol{t}_1/\gamma_1, \ldots, \boldsymbol{t}_m/\gamma_r)
$$

とすると，$(\boldsymbol{t}_i/\gamma_i)' (\boldsymbol{t}_j/\gamma_j)$ は $i=j$ のとき 1, さもなくば 0 である．

後半の例ではこれらの列（タテ）ベクトルを特異ベクトルと呼んでいる．

### 線形代数独学者のためのノート

上述した [佐竹『線形代数学』](https://www.shokabo.co.jp/mybooks/ISBN978-4-7853-1316-6.htm) では定理 2 と同じことが単に定理 1 の「例2」として述べられている（『Ⅳ 行列の標準化』の3節『対称行列の標準化』）．

線形代数は微分積分と並んで大学の理工系の学部の1〜2年生向けカリキュラムとされることが多い．つまり専門分野がはっきり決まる前の基本的な事柄というニュアンスが強いので，「特異値分解」という特別な用語を使わずもう少し一般的な形で記載されることが多いのだと思う．

とはいえ，なんの役に立つか（どういう応用があるか）わからないことを学ぶのは難しいので，一例として特異値分解みたいな応用があることを知っておくのは悪くないと思う．


## 数値計算してみる

R では特異値分解を行う関数が `svd` として提供されている．

`iris` というアヤメのガクと花弁の長さと幅を記録したデータを特異値分解してみる．

```r
X = as.matrix(iris[,-5])
muhat = colMeans(X)
X = sweep(X, 2, muhat)
res_svd = svd(X)
col3 =  hcl.colors(3)
plot(res_svd$u, col =col3[iris$Species], pch=16)
legend("topright", legend=levels(iris$Species), col=col3, pch=16)
```

![](/images/diagonal_to_svd/scatter_u.png)

ここでは特異値を大きい順に並べたときの上から2つめまでに対応する特異ベクトルをプロットした．

このようにデータをより小さい次元に縮約して図示することは，データの特徴を把握する目的でよく行われる．

ここでは例えば，versicolor と virginica は似ている個体もあるが，setosa と versicolor, virginica はあまり似ていないことがわかる．

次に，定理 2 を示したときのように，固有値・固有ベクトルから特異値分解を求めてみたい．
固有値・固有ベクトルを求める関数は `eigen` で提供されている（内部ではLAPACK という FORTRAN のプログラムを使っているそう）．


```r
res_eigen = eigen(t(X)%*%X)
```

ところで，行列 $A$ の固有値・固有ベクトルは

$$
A\boldsymbol{x} = \alpha \boldsymbol{x}
$$

を満たすものであった（再掲）．この式は両辺に同じ値を掛けても成り立つので固有値・固有ベクトルは一意でない．しかし，`eigen` では固有ベクトルの長さが 1 になるように正規化した値を出力するように約束してくれている．

そして異なる固有値に属する固有ベクトルは一次独立なので `eigen` の返してくれる固有ベクトル（`vectors`）は直交行列である．

```r
> round(t(res_eigen$vectors)%*%res_eigen$vectors, 5)
     [,1] [,2] [,3] [,4]
[1,]    1    0    0    0
[2,]    0    1    0    0
[3,]    0    0    1    0
[4,]    0    0    0    1
```

なので `vectors` をそのまま定理 2 の $P_2$ として使えば良い．

```r
P1 = sweep(X%*%res_eigen$vectors, 2, sqrt(res_eigen$values), FUN = "/")
```

確認すると次のように `svd` 関数と同じ値が求まっていることがわかる．


```r
 > print(res_svd$d)
 [1] 25.099960  6.013147  3.413681  1.884524
 > print(sqrt(res_eigen$values))
 [1] 25.099960  6.013147  3.413681  1.884524

 > print(res_svd$v)
            [,1]        [,2]        [,3]       [,4]
[1,]  0.36138659 -0.65658877  0.58202985  0.3154872
[2,] -0.08452251 -0.73016143 -0.59791083 -0.3197231
[3,]  0.85667061  0.17337266 -0.07623608 -0.4798390
[4,]  0.35828920  0.07548102 -0.54583143  0.7536574
> print(res_eigen$vectors)
            [,1]        [,2]        [,3]       [,4]
[1,]  0.36138659 -0.65658877  0.58202985  0.3154872
[2,] -0.08452251 -0.73016143 -0.59791083 -0.3197231
[3,]  0.85667061  0.17337266 -0.07623608 -0.4798390
[4,]  0.35828920  0.07548102 -0.54583143  0.7536574

```

```r
plot(P1, res_svd$u)
abline(0, 1, col="royalblue")
```

![](/images/diagonal_to_svd/comparison_u.png)


統計学の分野で用いられる主成分分析も特異値分解とほぼ同じものである．

定理 2 と同じ記号で書くと，主成分分析の場合は $(A P_2)'$ を主成分と呼ぶ．

主成分分析を行う関数は `prcomp` として提供されている．

```r
res_pca = prcomp(as.matrix(iris[,-5]))
```

特異値分解の特異値（`d`）をふたたび $P_1$（`u`）に掛けてやると主成分分析の結果と一致することがわかる．

```r
res_svd = svd(X)
PC = sweep(res_svd$u, 2, res_svd$d, FUN="*")

plot(PC, res_pca$x)
abline(0, 1, col="royalblue")
```

![](/images/diagonal_to_svd/comparison_pc.png)

他にも触れておきたい（触れておくべき？）ことはいっぱいあるが，それらについてはもう少し勉強してからでないと書けなかったりするし，ちょっと疲れてきたのでここで終わりにする．
