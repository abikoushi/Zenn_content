---
title: "TF-IDF指標とカルバック・ライブラ情報量"
emoji: "💭"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R]
published: false
---

## はしがき

Akiko Aizawa (2003), An information-theoretic perspective of tf–idf measures, Information Processing & Management,
(39) 1, https://doi.org/10.1016/S0306-4573(02)00021-3 

（の一部）を読んだメモです．ただし表記などは私の好みに合わせて変えてあります．自然言語処理などで使うTF-IDF指標がKLダイバージェンスの各項から出てくることなどが書かれています．

## 本文

$w$ と $d$ を有限の台を持つ離散型の確率変数とする．つまり多項分布を考える．$w$ は単語（ワード），$d$ は文書（ドキュメント）に対応する質的確率変数だと思ってほしい．

2つの分布 $f$, $g$ のカルバック・ライブラ（KL）情報量を $D[f \| g]$ と表記する．次の $\mathscr{M}$ は同時分布 $p(x,y)$ を独立な分布の積 $p(x) p(y)$ で近似したときのKL情報量である．

$$
\begin{aligned}
\mathscr{M} &= D[ p(w, d) ,\ p(w) p(d)) ] \\
&= \sum_{w,d} p(w,d) \cdot \log \frac{p(w,d)}{p(w) p(d) }.
\end{aligned}
$$

この $\mathscr{M}$ を期待相互情報量と呼ぶことにする．期待相互情報量は次のようにも書ける．

$$
\mathscr{M} = \sum_{w,d} \, p(w, d) \, \log \frac{p(d|w)}{p(d)}
\tag{1}
$$


さて, $f_{wd}$ を文書 $d$ における単語 $d$ の頻度の観測値とし，さらに $f = \sum _{w,d} f_{wd}$ , $f_{w} = \sum _d f_{wd}$ と書くことにする．

また， $D$ を文書の数（$d = 1, \ldots , D$），$D_w$ を単語 $w$ を含む文書の数とする．

この記号の下で，文書，単語をランダムに（つまり独立同分布で）選ぶとして， $p(w,d)$ の最尤推定量は，

$$
\hat {p}(w,d) = f_{wd}/f.
$$

$p(d|w)$ , $p(d)$ の最尤推定量はそれぞれ

$$
\hat{p}(d|w) = 1/D_w, \quad \hat{p}(d) = 1/D
$$

である．

最尤推定量を期待相互情報量 (1) の各項にプラグインすると，


$$
\hat {p}(w,d) \cdot \log \frac{\hat{p}(d|w)}{\hat{p}(d) }= \underbrace{\frac{f_{wd}}{f}}_{\text{TF}} \cdot \underbrace{\log \frac{D}{D_w}}_{\text{IDF}} 
$$

となる．右辺は TF-IDF（term frequency-inverse document frequency）と呼ばれ，文書ごとの特徴的な単語を抜き出す指標などとして用いられる．

以下は私見だが，この程度に抽象化しておくとちょっとした応用・変更が気軽にできる場合がある．例えば…

- ゼロ除算を避けるために IDF の分子分母に小さい数（1 など）を足すことは，最尤法でなく MAP 法を使うことに相当する
- 文書ごとでなく適当なクラスタリングで得られたクラスタごとにTF-IDF指標を求めることは，矛盾なくできる

最後に R による計算例を載せておく．

### 計算例

`janeaustenr` のデータを用いて TF-IDF 指標を計算してみる．Jane Austen 氏の著作ごとに特徴的な語を選ぶ．

```R
library(janeaustenr)
library(tidytext)
library(dplyr)
library(ggplot2)

#頻度ｆの集計
book_words <- austen_books() %>%
  unnest_tokens(word, text) %>%
  count(book, word, sort = TRUE) %>% 
  group_by(book) %>% 
  mutate(total = sum(n)) %>%
  ungroup() %>% 
  bind_tf_idf(term = word, document = book, n = n) %>% 
  arrange(desc(tf_idf))

#定義通りTF-IDFを求める
book_tf_idf = mutate(book_words, tf2 = n/total) %>% 
  group_by(word) %>% 
  mutate(idf2 = -log(n_distinct(book))) %>% 
  ungroup() %>% 
  mutate(idf2 = idf2 + log(n_distinct(book))) %>% 
  mutate(tf_idf2 = idf2*tf2)

##　tidytext パッケージの計算結果と一致することを確認
ggplot(book_tf_idf, aes(tf_idf2 , tf_idf))+
  geom_abline(slope = 1, intercept = 0, linetype=2)+
  geom_point(alpha=0.2)+
  theme_classic(14)

#大きい方から上位10個抜き出し
book_top <- book_tf_idf %>%
  group_by(book) %>%
  slice_max(tf_idf, n = 10) %>%
  ungroup()

#プロット
ggplot(book_top, aes(x = tf_idf, y = reorder_within(word, tf_idf, book))) +
  geom_col() +
  facet_wrap(~book, ncol = 2, scales = "free") +
  scale_y_reordered()+
  labs(x = "TF-IDF", y ="word") +
  theme_classic(14)
```

![](/images/tfidf_and_kld/tfidf_top10.png)

おしまい．