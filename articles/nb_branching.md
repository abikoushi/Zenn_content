---
title: "Rによる離散時間の分岐過程のシミュレーションとその感染症疫学への応用"
emoji: "🪾"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 疫学, 確率過程]
published: true
---

## このノートについて

このノートは [西浦 博; 小林 鉄郎; 安齋 麻美; 合原 一幸; ナタリー・リントン. 感染症流行を読み解く数理（日本評論社）](https://www.nippyo.co.jp/shop/book/8827.html) 第5章の「MERSの伝播リスク」の節の行間を補うためのものです．

同書の (5.9) 式のシミュレーションによる確認を行います．最初は導出も紹介しようと思ったのですが思ったより大変だったので挫折しました．導出については，

Nishiura, H., Yan, P., Sleeman, C. K., & Mode, C. J. (2012). Estimating the transmission potential of supercritical processes based on the final size distribution of minor outbreaks. Journal of theoretical biology, 294, 48-55. https://www.sciencedirect.com/science/article/pii/S0022519311005625?via%3Dihub

を見てください．

シミュレーションのコードは R です．

## シミュレーションしてみる

1人の感染者が次の感染者を生み出し，その感染者がまた次の感染者を生み出し，……というような過程を考えたい．

最初の1人の感染者を第0世代の感染者と呼ぶことにしよう．第0世代の感染者から生み出された感染者を第1世代，第1世代の感染者から生み出された感染者を第2世代，……，第 $n-1$ 世代の感染者から生み出された感染者を第 $n$ 世代の感染者と呼ぶことにする．

第 $n$ 世代の $i$ 番目の感染者が生み出す2次感染者数を $X_i(n)$ と置く．1人の感染者から生み出される感染者数 $X_i(n)$ を独立な負の二項分布とする．平均が $R_0$ となるよう，ここでの負の二項分布の確率関数は次のようにパラメータ化される．

$$
\Pr(X_i(n)=x) = \frac{\Gamma(k+x)}{x!\Gamma(k)}\left(\frac{R_0}{R_0+k}\right)^x \left(1+\frac{R_0}{k}\right)^{-k} \tag{1}
$$

ここでの $\Gamma(x)$ はガンマ関数である．この分布の平均は $R_0$，分散は $R_0+R_0^2/k$ となる

R の `*nbinom` 関数を使う際は，次のように置くと (1) の確率関数と一致する．

$$
\texttt{size} \leftarrow k, \quad \texttt{mu} \leftarrow R_0
$$

あるいは　$\texttt{prob} \leftarrow k/(k+R_0)$ としてもいい．

イメージをつかむため，まず9回のシミュレーション結果を図示してみる．9という中途半端な数なのは単に3×3のプロットに配置したかったからという理由による．$R_0=0.6$, $k=5$ とした．

![](/images/nb_branching/tree.png)
*横軸の感染者の添字（通し番号），縦軸が世代．なにもプロットされていないパネルは第1世代の感染者数が0だったことを示す．*

手元に『感染症流行を読み解く数理』をお持ちの方は図5.3と比較してみてほしい．図5.3のようなモデルがシミュレーションされていることがわかると思う．

シミュレーションを行う関数はモデルの仮定から愚直に次のように書いた．

```r
branching_process <- function(R0, k, G){
  gen <- vector("list", G)
  #p <- k/(k+R0)
  X <- rnbinom(1, size=k, mu=R0) #第1世代
  if(X>0){
    id <- 1L:X
    parent <- 0L
    gen[[1]] <- data.frame(id = id,
                           parent = parent,
                           g = 1L)
    parent <- id
    next_id <- X
    
    for(i in 2L:G){  #第2世代以降
      X <- rnbinom(length(parent), size = k, mu = R0)
      sumX <- sum(X)
      if(sumX == 0){ break } #感染者が0になったら終了
        parent <- rep(parent, X)
        id <- (next_id+1L):(next_id+sumX)
        gen[[i]] <- data.frame(id = id,
                               parent = parent,
                               g = i)
        parent <- id
        next_id <- next_id+sumX
    }    
  }
  dplyr::bind_rows(gen)
}
```

さて，結論のみ引用するが，第 $n$ 世代までの総感染者数の $n \to \infty$ の極限， $Y = \sum_{n=0}^{\infty}\sum_{i=1}X_i(n)$ の確率関数は次のようになる．

$$
\Pr(Y=y)=\begin{cases}
\frac{1}{(1+R_0/k)^k} & y=1\\
\frac{\prod_{j=0}^{y-2}(j/k+y)}{y!}\left(\frac{k}{R_0+k}\right)^{ky} \left(\frac{R_0k}{R_0+k}\right)^{y-1} & y\ge2
\end{cases}
$$

このことを10000回のシミュレーションで確認しよう．この総感染者数 $Y$ の確率関数は次のように実装した．

```r
prob_totalsize <- function(y,R0,k){
  if(y==1){
   1/(1+R0/k)^k
  }else{
    prod(0:(y-2L)/k+y)*((k/(R0+k))^(k*y))*(R0*k/(R0+k))^(y-1)/factorial(y) 
  }
}
```

無限世代までのシミュレーションは無理なので，ひとまず20世代まで行った．結果的に今回のシミュレーションでは最大15世代までで収束していることを事後的に確認した．

シミュレーションによる総感染者数の分布を次の図に示す．

![](/images/nb_branching/totalsize.png)
*縦棒がシミュレーションで得られた実現値，丸いマーカーが `prob_totalsize` 関数による理論値．*

シミュレーション結果が理論値と一致することがわかる．

作図も含めたRのコード全体はこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/nb_branching.R
