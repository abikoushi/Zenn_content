---
title: "フィッシャーの正確確率検定の紹介"
emoji: "🌊"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 仮説検定, 信頼区間]
published: false
---

## はじめに

フィッシャーの正確確率検定の紹介はWeb上にもいろいろなページがある. 

しかし [R の exact test についてのノート](https://zenn.dev/abe2/articles/exact_tests_r) のように

「帰無仮説のパラメータを引数にして p 値を返す関数が書かれていれば, こんなふうにして信頼区間も得ることができ, 2つが整合する」

という方針だと, 帰無仮説がフィッシャーの非心超幾何分布のパラメータになっていることを明示していない解説は「こっち見といて」としにくい. 

そこで自分でも書くことにした. 人によって強調するポイントが違ったりするのも悪くはないだろう.


## 超幾何分布

なんかの要因への曝露（ばくろ; 悪いこととは限らない）によって, ある疾病を持つかどうか調べたところ, 次のようなデータが得られたとしよう.

|| 曝露あり | 曝露なし |
| ----| ---- | ---- |
|病気あり|$x$ | $y$|
|病気なし|$m-x$ |$n-y$ |



この表が得られる確率を計算する一つの方法として,次のような2項分布のモデルが考えられる.

$$
x \sim \mathrm{Binom}(m, q)
$$

$$
y \sim \mathrm{Binom}(n, r)
$$

$k=x+y$ と置くことにし, 表の周辺度数（タテ・ヨコ合計）を固定したとき,

$$
\begin{aligned}
p(x|m,n,k,q,r) &\propto  p(k-x|n,r,x)p(x|m,q)\\
&\propto r^{k-x}(1-r)^{n-(k-x)} q^x(1-q)^{m-x}. \tag{1}
\end{aligned}
$$

$x$, $y$ それぞれのオッズを $\omega _x = q/(1-q)$, $\omega _y = r/(1-r)$, さらにオッズ比を $\omega=\omega_x/\omega_y$ と書くことにすると,

$$
\begin{aligned}
(1) &= \omega_y^{k-x}(1-r)^{n} \omega_x^x(1-q)^{m}\\
&=\omega^x \omega_y^{k}(1-r)^{n}(1-q)^{m}.
\end{aligned}
$$

また $k$ の分布は次のたたみ込みである.

$$
p(k|m,n,q,r) = \sum_u p(y=k-u|n,r,x)p(x=u|m,q).
$$

よって $x$ に依存しない $\omega_y^{k}(1-r)^{n}(1-q)^{m}$ の部分は分子・分母で消え,

$$
\begin{aligned}
p(x|m,n,k,q,r) &= \frac{p(k-x|n,r,x)p(x|m,q)}{p(k|m,n,q,r)}\\
&= \frac{1}{\sum_u \binom {m}{u} \binom {n}{n-u}\omega ^{u}} {\binom {m}{x}}{\binom {n}{n-x}}\omega ^{x}
\end{aligned}
$$

となる. これをフィッシャーの非心超幾何分布と呼ぶ（[Fisher's noncentral hypergeometric distribution - Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_noncentral_hypergeometric_distribution)）.

右辺を見ると $q$, $r$ を指定しなくてもオッズ比 $\omega$ を決めれば分布が一意に決まることがわかる.

上では $\sum_u$ の $u$ の範囲をあいまいに書いてしまったが, 注意が必要な場合もあるので注意してほしい.

また, 一般に確率の比が一致する確率分布は一致するの二項分布・二項分布でなく分割表のセルごとに独立なポアソン分布やカテゴリ数が $4=2\times2$ の多項分布で考えても同じになる.

##

(Polya's urn の特殊な場合).

|| 黒玉 | 白玉 |
| ----| ---- | ---- |
|取り出した|$x$ | $k-x$|
|取り出してない|$m-x$ |$n-(k-x)$ |


周辺度数（表の縦・横合計）を固定すると, 取りうる値の範囲をすべて書き出すことができ, 表が得られる確率を計算できる.

![](/images/hypergeometric/note_hypergeo.jpg)

（画像は「全部書き出して計算」の例）

この確率分布を超幾何分布という.

左上のセルの $x$ に注目すると周辺度数で条件付けることで $2 \times 2$ の分割表が1次元で表現できていることがわかる.


## フィッシャーの非心超幾何分布

超幾何分布は黒玉も白玉も平等な確率で取り出すとしたが, どちらかがより取り出しやいような場合も含むよう拡張したのがフィッシャーの非心超幾何分布（[Fisher's noncentral hypergeometric distribution - Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_noncentral_hypergeometric_distribution)）である. 通例, つぎのようにオッズ比 $\omega$ を使ってパラメトライズする.

$$
p(x) =  {\frac {{\binom {m_{1}}{x}}{\binom {m_{2}}{n-x}}\omega ^{x}}{P_{0}}}
$$

ここで $P_{0}$ は規格化定数（ここでは省略）.

この導出については [Fisher's noncentral hypergeometric distribution - Wikipedia](https://en.wikipedia.org/wiki/Fisher%27s_noncentral_hypergeometric_distribution) にも書かれているので, 
ここではちょっと楽をして黒玉・白玉のシミュレーションがこのフィッシャーの非心超幾何分布になることをたしかめてみよう.

（個人的には）
