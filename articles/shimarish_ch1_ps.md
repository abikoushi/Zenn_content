---
title: "しまりす本のP値関数とS値関数をRで再現する"
emoji: "🐿️"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R]
published: true
---

## 本文

タイトルで「しまりす本」と呼んでいるのは [佐藤俊哉『宇宙怪人しまりす統計より重要なことを学ぶ』（朝倉書店）](https://www.asakura.co.jp/detail.php?book_code=12297&srsltid=AfmBOooOBJ6JBrHb5K8-c_3nrza42gfLwEjXY5knjbjFV6ju-KR_tHHT) です．以下でもしまりす本と略記します．

まだ途中までしか読めてないけど，かなり良さそうです．

しまりす本の紹介は別の機会にあらためて書きたいのですが，本記事はタイトルの通り，単にこの本の第1話で出てくるP値関数とS値関数の図をR言語で再現してみるというものです．

分析対象は下の表のような，ヨクナール（架空の薬）の使用者・未使用者に対して5日目時点で風邪から回復した人の数を調べたデータです．

||回復|未回復|
|--|--|--|
|ヨクナール使用|58|22|
|経過観察（ヨクナールなし）|62|38|

しまりす本の書き方は数学的になりすぎないように気を使ったスタイルだと思いますが，このノートでするのは細かい話だけなので少し記号を使います．先の表に次のように記号を対応させます．

||回復|未回復|
|--|--|--|
|$i=1$|$y_1$|$n_1 - y_1$|
|$i=2$|$y_2$|$n_2 - y_2$|

「ヨクナール使用」（$i=1$）と「経過観察」（$i=2$）の回復割合の推定量をそれぞれ $\hat{p}_i = y_i/n_i$ ($i=1,2$) とします．

推定量の標準偏差を標準誤差と呼びます．しまりす本では，「標準誤差は疫学や統計学の教科書に出ているから、それを使って計算すると、」とさらっと述べられていますが，複数のやり方があるので少し補足します．

$y_i$ が試行回数 $n_i$, 成功確率 $p_i$ の二項分布に従うとすると，二項分布の性質から，$\hat{p}_i$ の分散 $V[p_i]$ は，

$$
\begin{aligned}
V[\hat{p}_i] &= V[y_i/n_i]\\
&=\frac{1}{n_i^2} V[y_i] \\
&=\frac{1}{n_i^2} n_i p_i(1-p_i)\\
&=\frac{1}{n_i}p_i(1-p_i)
\end{aligned}
$$

となります．

回復割合の差 $\beta = p_1 - p_2$ についての推定量 $\hat \beta$ の分散は，分散の加法性から，

$$
\begin{aligned}
V[\hat \beta] &= \sum_{i=1}^2 V[y_i/n_i]\\
&= \sum_{i=1}^2\frac{1}{n_i}p_i(1-p_i)
\end{aligned}
$$

です．よって標準誤差 $se[\hat \beta]$ は分散の平方根を取って，

$$
se[\hat \beta] = \sqrt{\sum_{i=1}^2\frac{1}{n_i}p_i(1-p_i)}
$$

です．が，$p_i$ は不明なので標準誤差についてもなんらかの推定をします．標準誤差そのものでなく標準誤差の推定量を使うので複数のやり方がありえます．

しまりす本で採用しているのは $p_i$ に $\hat p_i$ をそのまま代入する方式です．

一方で，教科書によっては　$p_i$　に次の $\hat p_{\text{pooled}}$ を代入する方式が紹介されている場合もあります．

$$
\hat p_{\text{pooled}} = \frac{\sum_{i=1}^2 y_i}{\sum_{i=1}^2 n_i}.
$$

これは $\beta=0$ とした場合，$p_{\text{pooled}} = p_1 = p_2$ なので，どうせなら全データ使って推定しようというものです．今回は差があるとした帰無仮説も含めて計算したいので1個目の標準誤差のほうが便利です．

ここからRによる実装例です．

先に「差がない」とした2個目の方の推定値で求めた標準誤差を記載します．これはしまりす本に記載されている標準誤差「0.0696」と少し違います．

```r
#データ
X = matrix(c(58, 22,
             62, 38), byrow = TRUE, nrow = 2, ncol = 2)
print(X)
####
#「差がない」とした推定値で求めた標準誤差
phat = x/n
p_pooled = sum(x)/sum(n)
delta = phat[1]-phat[2]
se0 = sqrt(p_pooled*(1-p_pooled)*sum(1/n))
print(se0)
#[1] 0.07071068
```

一方でP値を見ると，「13.8％」という記述と一致しています．

```r
#自由度1のカイ二乗分布
print(pchisq(abs(delta/se0)^2, 1, lower.tail = FALSE))
#[1] 0.1375639
```

これもこのノートを書こうと思ったきっかけの一つで，要はしまりす本の「13.8％」というのが誤植で次の「13.2％」が正しいと思われます（私が持ってるのはKindle版の初版）．

[経緯（リンク先Xの投稿）:黒木玄 Gen Kuroki (@genkuroki), July 8, 2025 ](https://twitter.com/genkuroki/status/1942441357248979484?)

1個目の差を許すの方の推定値で標準誤差を求めると「0.0696」と一致します．

```r
#差を許す推定値で求めた標準誤差
se = sqrt(sum(phat*(1-phat)/n))
print(se)
#[1] 0.06962893

pvalfun_wald <- function(mu0, delta, se){
  pchisq(((delta-mu0)/se)^2, df = 1, lower.tail = FALSE)  
}

cat("H0:diff=0",round(pvalfun_wald(0, delta=delta, se=se)*100, 1),"%")
#H0:diff=0 13.2 %
```

さて，準備が整ったのでP値関数をプロットします．

```r
pvalfun_wald <- function(mu0, delta, se){
  pchisq(((delta-mu0)/se)^2, df = 1, lower.tail = FALSE)  
}

mu_v <- seq(-0.2, 0.4, length.out=501)
res_p = sapply(mu_v, pvalfun_wald, delta=delta, se=se)

#png("Pfun.png")
plot(mu_v, res_p, type="l", xlab = "risk difference", ylab = "p-value")
#dev.off()
```

![](/images/shimarish_ch1_ps/Pfun.png)

これがやりたかったことでした．

これで完成でもいいのですが，信頼区間も書き込んで，もうちょっと本に出てるような図に寄せてみます．

```r
library(ggplot2)
df_p = data.frame(hypo=mu_v, pv=res_p)

cifun_wald <- function(level, delta, se){
  alpha = (1-level)*0.5
  z = qnorm(1-alpha,0,1)
  data.frame(level=level, lower=delta - se*z,   upper=delta + se*z)
}

ci1 = cifun_wald(c(0.8,0.95),delta,se)

p1 = ggplot(df_p)+
  geom_line(aes(x=hypo, y=pv))+
  geom_errorbarh(data = ci1, 
                 aes(xmin = lower, xmax=upper, y=1-level, colour=factor(level)), 
                 height = 0.05)+
  geom_vline(xintercept = 0, colour="grey")+
  scale_color_brewer(palette = "Set1")+
  theme_bw(18, base_family = "Osaka") + 
  labs(x="設定した回復割合の差", y="両側P値", colour="信頼係数")
print(p1)

#ggsave(filename = "Pfun_gg.png", plot = p1)
```

![](/images/shimarish_ch1_ps/Pfun_gg.png)


私の考えとしては，

- 「仮説検定って有意水準ちょっと変えたらたらどうなるの？　どうせなら全部やろう」がP値
- 「仮説検定って帰無仮説ちょっと変えたらどうなるの？　どうせなら全部やろう」が信頼区間

と捉えると，P値関数（＝信頼区間関数）を見たい気持ちが自然に思えてきます．どっちも全部やったものがP値関数だからです．

S値はP値にたいして底が 0.5 の対数を取ったものです．せっかくなので同様にプロットしておきます．

```r
df_s = data.frame(hypo=mu_v, sv=log(res_p,base = 0.5))
p2 = ggplot(df_s)+
  geom_line(aes(x=hypo, y=sv))+
  geom_vline(xintercept = 0, colour="grey")+
  theme_bw(18, base_family = "Osaka") + 
  labs(x="設定した回復割合の差", y="ビット")
print(p2)

#ggsave(filename = "Sfun_gg.png", plot = p2)
```

![](/images/shimarish_ch1_ps/Sfun_gg.png)

できました．




## 落穂拾い

ちなみに [ロジスティック回帰を経由して2×2の分割表のオッズ比の検定を作る](https://zenn.dev/abe2/articles/5ef89a9f5b2ab6) のような形でモデルを書くと「回復割合の差」は次の $\beta$ に相当し，恒等リンク関数を用いた回帰モデルで考えるのと同等とわかる．

$$
\begin{aligned}
y_i &\sim \mathrm{Binom}(n_i, p_i)\\
& \text{where}\quad p_i = \alpha + \beta x_i\\
& s.t. \quad 0\le p_i \le 1 \quad (i= 1,2)
\end{aligned}
$$


また，上記で紹介したような $T(X) = |\hat \theta - \theta| / se(\hat \theta)$ で検定統計量を与えるタイプの検定をワルド（Wald）検定と呼ぶ（[尤度比検定, ワルド検定, スコア検定から定まる信頼区間：二項分布の例](https://zenn.dev/abe2/articles/lr_wald_score_binom) ）．

割合の差のスコア検定については，

- [割合の差のスコア検定 - Triad sou.](https://triadsou.hatenablog.com/entry/2024/03/01/001513)
- [母比率の差に関するP値と信頼区間のスコア法による構成 - 黒木玄](https://nbviewer.org/github/genkuroki/public/blob/main/0047/score%20method%20for%20risk%20difference.ipynb) 

に解説がある．

