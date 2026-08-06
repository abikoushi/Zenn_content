---
title: "Rによる連続時間分枝過程のシミュレーション"
emoji: "🚸"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 確率過程]
published: true
---

## はじめに

[近江崇宏・野村俊一『点過程の時系列解析』（共立出版）](https://www.kyoritsu-pub.co.jp/book/b10003181.html) の5章4.2節で，次のような時間構造を取り入れた分枝過程（branching process; 分岐過程）が紹介されている．

1. 第1世代のノードは強度 $\mu$ の定常ポアソン過程に従って発生する．
2. 時刻 $t_i$ に生成されたノードは，強度関数 $g(t-t_i)$ ($t>t_i$) の非定常ポアソン過程に従って次の世代の子ノードを生成する．

この分枝過程のシミュレーションをやってみたいが，今回は少し設定を単純化し，次の確率過程を考えることにする．

1. 第1世代のノードは強度 $\mu$ の定常ポアソン過程に従って発生する
2. 各ノードは，強度 $\lambda$ の定常ポアソン過程に従って次の世代の子ノードを生成する

この場合，ある時刻 $t$ までに生成されたノードの数の期待値は陽に書き表わすことができる．

$Z_n(t)$ を時刻 $t$ までに出生した第 $n$ 世代の個体数とする．全個体数 $N(t)$ は次のように表せる．

$$
N(t)=\sum_{n=1}^{\infty}Z_n(t)
$$

まず，第1世代のノードの数はパラメータ $\mu t$ のポアソン分布に従う．

$$
Z_1(t) \sim \mathrm{Poisson}(\mu t).
$$

第1世代のノードの数の期待値は $m_1(t) = E[Z_1(t)]=\mu t$ である．

第2世代のノードの数の期待値 $m_2(t)$ は，

$$
m_2(t)= 
\int_0^t \lambda  m_1(s) \,ds =\frac{\mu\lambda t^2}{2}.
$$

第3世代のノードの数の期待値 $m_3(t)$ は，

$$
m_3(t)=\int_0^t \lambda m_2(s)\,ds =
\frac{\mu \lambda^2 t^3}{3!}.
$$

第 $n$ 世代のノードの数の期待値は 

$$
m_n(t)=\frac{\mu\lambda^{n-1}t^n}{n!}.
$$

したがって全個体数の期待値は世代ごとの期待値を全世代にわたって合計し，次のようになる．

$$
\begin{aligned}
E[N(t)]&= \sum_{n=1}^{\infty}
\frac{\mu\lambda^{n-1}t^n}{n!}\\
&= \frac{\mu}{\lambda}\sum_{n=1}^{\infty}
\frac{\lambda^n t^n}{n!}\\
&=\frac{\mu}{\lambda}
\left(
e^{\lambda t}-1
\right)
\end{aligned}
$$

最後の等号は指数関数のテイラー展開による．

## シミュレーション

シミュレーションのコードは上の2つの仮定から愚直に，次のように書いてみた．

```r
simulate_branching <- function(mu,
                               lambda,
                               Tmax) {
  N0 <- rpois(1, Tmax*mu) #第1世代のノードの数
  t <- sort(runif(N0, 0, Tmax)) #発生時刻をポアソン配置
  individuals <- data.frame(
    id = 1L:N0,
    parent = 0L,
    birth_time = t,
    generation = 1L
  )
  
  next_id <- N0 + 1L
  i <- 1L
  while (i <= nrow(individuals)) {
    birth <- individuals$birth_time[i]
    t <- birth
    repeat {
      t <- t + rexp(1, rate = lambda) #子ノードの発生時刻
      if (t > Tmax){
        break
      }
      #Tmax 未満ならノードを追加
      individuals <- rbind(
        individuals,
        data.frame(
          id = next_id,
          parent = individuals$id[i],
          birth_time = t,
          generation = individuals$generation[i] + 1L
        )
      )
      next_id <- next_id + 1L
    }
    i <- i + 1L
  }
  
  individuals <- individuals[order(individuals$birth_time), ]
  return(individuals)
}
```

R ではこのように行を逐次的に追加していくコードは遅いことがよく言われるが，今回はこれでいいことにする．

イメージをつかむため，まず1回のシミュレーション結果を図示してみる．この図は『点過程の時系列解析』の図5.4を参考にしたものである．手元に『点過程の時系列解析』をお持ちの方は図5.4と比較してみてほしい．

![](/images/sim_branchingp/tree.png)
*分枝過程のシミュレーション結果．　$\mu=1$, $\lambda=0.6$とした．*

次に時間ごとのノードの数を示す．これは100回のシミュレーションの結果をまとめて図示する．

![](/images/sim_branchingp/count.png)
*分枝過程の時間ごとのノードの数．青い点線は平均．*

ポアソン過程と比較すると，ばらつきが大きいことがわかる．


作図も含めたコード全体は以下に置く：

https://github.com/abikoushi/Zenn_content/blob/main/R/sim_branchingp.R