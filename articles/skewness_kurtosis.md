---
title: "Julia で学ぶ歪度・尖度"
emoji: "🎃"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [Julia, 確率分布]
published: false
---

## あらまし

モーメント母関数の対数を取ったものをキュムラント母関数という．
キュムラント母関数のテイラー展開を用いると中心極限定理による正規分布への収束の速さが歪度・尖度でほぼ決まることがわかる．

また，Julia の [Distibutions.jl](https://juliastats.org/Distributions.jl/stable/) パッケージでは確率分布の歪度・尖度の計算も実装されているので，これを使って数値的にもたしかめてみる．

## 証明してみる 

このパートは [竹内啓『数理統計学: データ解析の方法』（東洋経済新報社）](https://books.google.co.jp/books/about/数理統計学.html?id=fSo8DwAAQBAJ&redir_esc=y) を参考にしている．

### 準備

確率変数 $X$ の **モーメント母関数** $\phi_X(t)$ を

$$
\phi_X(t) = E[e^{tX}]
$$

と定義すると，テイラー展開したときに


### 

## 計算してみる

