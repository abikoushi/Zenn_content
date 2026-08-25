---
title: "分枝過程と Hawkes 過程の関係"
emoji: "🥦"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 確率過程, 統計学]
published: true
---

## はじめに

分枝過程（branching process; 分岐過程）は Hawkes 過程より解析的な扱いが容易であり，分枝過程を解析することで Hawkes 過程の性質を調べることがあるそうだ．

[近江崇宏・野村俊一『点過程の時系列解析』（共立出版）](https://www.kyoritsu-pub.co.jp/book/b10003181.html) の5章4.2節で，時間構造のある分枝過程から発生したノードの発生時刻は Hawkes 過程に従うことが述べられている．

数学的な証明は『点過程の時系列解析』にあたってもらうことにして，ここでは R 言語による数値例を通じてこのことを見ていく．

区間 $[0, T]$ で次の確率過程を考える．

1. 第1世代のノードは強度 $\mu$ の定常ポアソン過程に従って発生する
2. 各ノードは，パラメータ $\beta$ の指数分布に従う生存時間を持つ
3. 各ノードは，生存時間の間だけ強度 $\alpha$ の定常ポアソン過程に従って次の世代の子ノードを生成する

このとき，ノードの発生時刻は指数カーネルの Hawkes 過程に従う．

ノード $i$ の寿命を $L_i$ とし，これがパラメータ $\beta$ の指数分布に従う．

$$
L_i \sim \operatorname{Exp}(\beta)
$$

発生時刻 $t_i$ のノード $i$ は $t_i+L_i$ で死亡し，それ以降は子ノードを産まなくなる．

$L_i$ を所与としたときのノード $i$ の子ノードの発生する強度は，

$$
g_i(t|L_i) =\alpha \mathbf 1(t_i<t<t_i+L_i)
$$

ここで $\mathbf 1(t_i<t<t_i+L_i)$ は指示関数とした（$t_i<t<t_i+L_i$ を満たすとき1，そうでなければ0の値を取る）．

指数分布の仮定から，

$$
P(L_i>t-t_i)=
e^{-\beta(t-t_i)}.
$$

よって，ノード $i$ の子ノードの発生する強度の期待値は，

$$
g_i(t) = E[g_i(t|L_i)]=\alpha e^{-\beta(t-t_i)}.
$$

これを過去の履歴 $H_t$ について足して，時刻 $t$ における条件付き強度は，

$$
\lambda(t|H_t)=
\mu+
\alpha
\sum_{t_i<t}
e^{-\beta(t-t_i)}
$$

となる．これは指数カーネルを用いた場合の Hawkes 過程の条件付き強度関数と一致する．

## とりあえずシミュレーションしてみる

上記の分枝過程をシミュレーションする関数を次のように書いた．

```r
simulate_branching_death <- function(mu,
                                     alpha,
                                     beta,
                                     Tmax) {
  
  N0 <- rpois(1, Tmax * mu)
  if(N0==0L){
    return(
      data.frame(
        id = NA,
        parent = 0L,
        birth_time = NA,
        death_time = NA,
        generation = NA
      )
    )
  }
  
  t0 <- sort(runif(N0,0,Tmax))
  individuals <- data.frame(
    id = seq_len(N0),
    parent = 0L,
    birth_time = t0,
    death_time = t0 + rexp(N0, rate = beta),
    generation = 1L
  )
  
  next_id <- N0 + 1L
  i <- 1L
  
  while (i <= nrow(individuals)) {
    
    birth <- individuals$birth_time[i]
    death <- min(individuals$death_time[i], Tmax)
    
    if (birth < death) {
      
      t <- birth
      
      repeat {
        
        t <- t + rexp(1, rate = alpha)
        
        if (t > death || t > Tmax) {
          break
        }
        
        individuals <- rbind(
          individuals,
          data.frame(
            id = next_id,
            parent = individuals$id[i],
            birth_time = t,
            death_time = t + rexp(1, rate = beta),
            generation =
              individuals$generation[i] + 1L
          )
        )
        
        next_id <- next_id + 1L
      }
    }
    
    i <- i + 1L
  }
  return(individuals[order(individuals$birth_time),])
}
```

$T = 10$, $mu = 0.5$, $a = 0.8$, $b = 1$ とパラメータを指定し，シミュレーション結果を図示してみよう．ちょっとごちゃごちゃした図になるがシミュレーションしたい確率過程をイメージしやすいよう，ノードの死亡までの時間も描き込んだ．

![](/images/hawkes_and_branching/birth_death.png)
*縦軸が世代，横軸が時間，矢印が親子関係，丸がノードの発生，バツがノードの死亡を表す．*

各ノードの親子関係の情報を無視してノードの個数だけを示してみよう．

![](/images/hawkes_and_branching/birth_death_marginal.png)
*Birth がその時点までに発生したノードの数，Currentがその時点にすで発生していてまだ死亡していないノードの数，Death がその時点までに発生したノードの数を表す*

このノードの発生時刻が Hawkes 過程に従うということである．

## 最尤推定

分枝過程から発生したノードの発生時刻は Hawkes 過程に従うということはすなわち，ノードの発生時刻だけが観測されている状況をHawkes 過程とみて，分枝過程のパラメータの最尤推定量を求めることができる．

パラメータをまとめて $\theta = (\mu, \alpha, \beta)'$，ノードの総数を $n$ として，対数尤度関数は次のように書ける．

$$
\begin{aligned}
\log L(\theta)=\sum_{i=0}^{n-1}&\log \left[ \mu+ \sum_{j<i}\alpha \exp(-\beta(t_i-t_j)) \right]\\
&-\left[ \mu T + \sum_{i=1}^n \alpha(1-\exp(-beta(T-t_i))) \right].
\end{aligned}
$$

この対数尤度関数を R で実装すると次のようになった．

```r
ll <- function(par, ti, Tau){
  mu <- exp(par[1])
  b <- exp(par[2])
  a <- exp(par[3])
  sum(log(mu+sapply(ti, function(t)sum(a*exp(-b*(t-ti[ti<t])))))) - 
    (mu*Tau+sum(-a*expm1(-b*(Tau-ti))/b))
}
```

最尤推定量はこの関数を最大化する `par` を求めることで得られる．その最適化は関数 `optim` に丸投げできる．

先のシミュレーション・データから最尤推定を行い，シミュレーションで設定した真値から得られた条件付き強度関数と，最尤推定されたパラメータから得られた条件付き強度関数を比較してみよう．

![](/images/hawkes_and_branching/intensity.png)
*上段が累積のノード数と累積強度関数，下段がノードの発生時刻と強度関数．*

最尤推定値を使用したほうが手元のデータにはフィットしていることがうかがえる．

シミュレーションを1000回繰り返して最尤推定量の分布をもとめてみる．次のヒストグラムでは推定値と真値の差をプロットしている．

![](/images/hawkes_and_branching/hist_mle.png)
*シミュレーション結果．上から $\mu$,$\alpha$,$\beta$ について，推定値-真値のヒストグラム．*

推定値と真値の差なので0に近いほど嬉しい．$\mu$ に関してはそこそこうまく求まっている．$\alpha$,$\beta$ に関してはけっこう極端な値がまじってしまうことがわかる．

## 分枝比

分枝比（branching ratio）と呼ばれる量 $\gamma$ も見ておこう．分枝比はカーネル関数 $g_i(t)$ を区間 $[0, \infty]$ で積分した値として定義される．

$$
\begin{aligned}
\gamma &= \int_0^\infty
g_i(u) \, du\\
&=\int_0^\infty\alpha e^{-\beta u} \, du\\
&=\frac{\alpha}{\beta}.
\end{aligned}
$$

これは感覚的にもわかりやすい結果である．なぜならば，平均寿命 $1/\beta$ のノードがその間に単位時間あたり $\alpha$ の子ノードを産むので，期待される1ノードあたりの子ノードの数

$$
E[\text{子ノードの数}]
= \alpha E[L]=\frac{\alpha}{\beta}
$$

とも一致しているからである．

分枝比は Hawkes 仮定および分枝過程を特徴づける量であり，

$$
\gamma<1
$$

ならばノードの数は平均的には $T\to\infty$ で有限の値に収束し，

$$
\gamma>1
$$

ならば平均的には $T\to\infty$ でノードが無限に増え続ける．

これまでのシミュレーションでは分枝比が0.8だったので分枝比が1より大きい場合もシミュレーションしてみよう．$\alpha=1.2$ とした．

次の図では生存時間の情報は省略する．

![](/images/hawkes_and_branching/birth_death2.png)
*縦軸が世代，横軸が時間，矢印が親子関係，丸がノードの発生を表す．*

ノードの個数は次のようになった．

![](/images/hawkes_and_branching/birth_death_marginal2.png)
*Birth がその時点までに発生したノードの数，Currentがその時点にすで発生していてまだ死亡していないノードの数，Death がその時点までに発生したノードの数を表す*

この記事を書くのに使用したコードのまとめはこちら：

https://github.com/abikoushi/Zenn_content/blob/main/R/hawkes_and_branching.R