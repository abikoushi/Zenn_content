---
title: "しまりす本第1話より：研究結果の再現性を分割表で考える"
emoji: "🌰"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計学]
published: false
---

## あらまし

タイトルの「しまりす本」というのは [宇宙怪人しまりす 佐藤俊哉『統計よりも重要なことを学ぶ』（朝倉書店）](https://www.asakura.co.jp/detail.php?book_code=12297&srsltid=AfmBOoqWTBxVz04vMTbTXhM3h4aslCpQeqNDTU7t9xscd3YlMcfydn-3) のことです．以下でも「しまりす本」と略します．

しまりす本第1話前半の要点は，

- 統計的仮説検定が完全に理想的に使われていたとしても，本当に効果がある候補が少ないときは再現性は低い
- 不正をしてる悪いやつがいるという発想ではこれは解決しない

です．このことを数値的に確認してみます．

## 準備

薬の候補となる物質のうち真に効果があるものを仮説検定で探すシチュエーションを考える．確率変数 $S$ を真に効果があるとき $S=1$， 真に効果がないとき $S=0$ の値を取るものとする．

また，確率変数 $T$ を仮説検定で有意になったとき $T=1$，有意でないとき $T=0$ の値を取るものとする．

仮説検定のαエラーというのはそっけない名前だが，真に効果がないときに検定で有意になる確率を表す．つまり，確率 $\alpha$ は次である．

$$
\alpha=P(T=1|S=0).
$$

仮説検定のβエラーというのもまたそっけない名前だが，真に効果があるときに検定で有意にならない確率を表す．つまり，確率 $\beta$ は次である．

$$
\beta=P(T=0|S=1).
$$

慣例的な仮説検定の枠組みの一つでは，αエラーを5％とか1％に固定した上で，検出力 $1-\beta$ が80％とかの目標に達するようサンプルサイズを決める．

このとき検定で有意になったものが真に効果がある確率は，次の式で表せる．

$$
\begin{aligned}
P(S=1|T=1) &= \frac{P(S=1, T=1)}{P(T=1)}\\
&= \frac{P(T=1 |S=1)P(S=1)}{P(T=1 ,S=1)+P(T=1 ,S=0)}\\
&= \frac{P(T=1 |S=1)P(S=1)}{P(T=1 |S=1)P(S=1)+P(T=1 |S=0)P(S=0)}\\
&=\frac{(1-\beta) P(S=1)}{(1-\beta) P(S=1)+\alpha P(S=0)}
\end{aligned}
$$

## Rによる数値例

上記 $P(S=1|T=1)$ の最後の式をそのまま打ち込む．

```r
reproducibility = function(true_ratio,
                           beta = 0.2,
                           alpha = 0.05){
  ((1-beta)*true_ratio)/((1-beta)*true_ratio+alpha*(1-true_ratio))
}
```

適当な値をいくつか与えてプロットしてみよう．

```r
library(dplyr)
library(ggplot2)
z = seq(0, 1, by=0.01) #P(S=1)
condf = expand.grid(alpha = c(0.01, 0.05, 0.1), 
                    beta = c(0.1,0.2,0.3)) %>% 
  group_by(alpha, beta) %>% 
  reframe(z=z,
          value = sapply(z, reproducibility, alpha=alpha, beta=beta))


ggplot(condf, aes(x=z,
                  y=value))+
  geom_line() + 
  facet_grid(alpha~beta, labeller = label_both)+
  labs(x = "true positive rate in the population",
       y = "true positive rate in the test positive")+
  scale_y_continuous(limits=c(0,1))+
  theme_bw(16)
```

![](/images/shimarisu_ch1_repro/TPR.png)
*横軸が $P(S=1)$，縦軸が $P(S=1|T=1)$*

αエラーとβエラーを一定の水準にたもっていたとしても，本当に効果のある候補が少ないときは $P(S=1|T=1)$ も低いことがわかる．

なので仮説検定で有意になったからといって正しい保証が得られたわけではない，背景の情報も加味して総合的に判断しましょうということになる．

P値を有意・有意でないの二分法でなく，連続的指標として使う方風については [しまりす本のP値関数とS値関数をRで再現する](https://zenn.dev/abe2/articles/shimarish_ch1_ps) も参照してほしい．

## おまけ：Rによる数値例2

この記事の「準備」でやったような条件付き確率の計算は地味にややこしく，私もよく間違える．そのため，答え合わせの方法として，次のようなモンテカルロ積分をやっておくと安心する．

```r
n_item <- 100000
p_population <- 0.05 #true positive rate in population
p <- c(0.05, 0.8) #alpha, 1-beta
set.seed(123)
pop = rbinom(n_item, 1, p_population)
test1 = rbinom(n_item, 1, p[pop+1])

true_positive = sum(test1 == 1 & pop == 1)/sum(test1 == 1)
```

割と仮定そのままを素直に計算できる感じがしないだろうか．

結果，下のように近い値が得られることがわかる．

```r
> print(true_positive)
[1] 0.4577594
> print(reproducibility(true_ratio = p_population, alpha=p[1], beta = 1-p[2]))
[1] 0.4571429
```

以上です．