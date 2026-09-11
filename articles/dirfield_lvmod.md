---
title: "R による方向場の描き方（ロトカ・ヴォルテラ方程式を例に）"
emoji: "🦈"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 微分方程式]
published: true
---

## 本文

次の連立微分方程式を考えます．

$$
\begin{aligned}
x' &= ax - bxy\\
y' &= -cy + dxy
\end{aligned}
$$

この微分方程式をロトカ・ヴォルテラ方程式と呼びます．$x$ は被食者（prey; 例えばイワシ），$y$ は捕食者（predator; 例えばサメ）の個体数だと思うとイメージしやすいです．各パラメータは次のような意味を持ちます．

- $a$ : $x$ の出生率
- $b$ : 被食による $x$ の死亡率
- $c$ : $y$ の死亡率
- $d$ : 捕食による $y$ の増加率

パラメータ $a$, $b$, $c$, $d$ はすべて正とします．

この微分方程式を R 言語で数値的に解くには次のようにします．

```r
library(dplyr)
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dPrey        <- a*Prey - b*Prey*Predator
    dPredator    <- -c*Predator + d*Prey*Predator
    return(list(c(dPrey, dPredator)))
  })
}

pars  <- c(a = 1, 
           b = 2, 
           c = 1,
           d = 3)

yini  <- c(Prey = 1, Predator = 0.5)
times <- seq(0, 50, by = 0.1)
out   <- ode(yini, times, LVmod, pars)
```

解をプロットしてみます．

![](/images/dirfield_lvmod/lvmod1.png)

タイトルにある方向場とは下の図のように微分方程式の解の方向を $(x,y)$ の平面上に示したものを指します．

![](/images/dirfield_lvmod/lvmod_dir.png)

図を見ると点線またぐところで方向が性質的に大きく変わることがわかります．左上から反時計回りに，

- $x$: 減少, $y$: 減少
- $x$: 増加, $y$: 減少
- $x$: 増加, $y$: 増加
- $x$: 減少, $y$: 増加

です.

点線は $(x', y')$ が 0 になる点の集合で,これをヌルクラインと呼びます．

今回のロトカ・ヴォルテラ方程式の場合は，$y=a/b$ と $x=c/d$ がヌルクラインです．

方向場の図は以下のような方法で描きました．

まず，いろいろな初期条件 `df_state` で微分方程式を片っ端から解く関数を宣言します．

```r
sol_ode_from_states <- function(df_state, times, func, parms){
  res_df <- lapply(1:nrow(df_state), function(i){
    ode_out <- ode(y=c(unlist(df_state[i,,drop=FALSE])), times=times, func=func, parms=parms)
    data.frame(group=i, ode_out)
  })
  dplyr::bind_rows(res_df)
}
```

初期値を指定して短時間だけ解きます．`expand.grid` は入力されたベクトルの要素のすべての可能な組み合わせを列挙したデータフレームを返す関数です．

```r
df_state <- expand.grid(Prey=seq(0.1, 2, 0.1),
                        Predator=seq(0.1, 2, 0.1))

times <- seq(0, 0.1, by=0.02)
res <- sol_ode_from_states(df_state=df_state, times=times, func=LVmod, parms=pars)
```

プロットはこちら．

```r
ggplot(data=res, aes(x=Prey, y=Predator, group=group, color=time)) +
  geom_path(arrow = arrow(length = unit(3, "pt"))) +
  scale_colour_gradient2(low="grey80", high="steelblue")+
  geom_hline(yintercept = pars["a"]/pars["b"], linetype=2)+
  geom_vline(xintercept = pars["c"]/pars["d"], linetype=2)+
  theme_bw()
```

R のコードを全部まとめたものはこちらです：
https://github.com/abikoushi/Zenn_content/blob/main/R/dirfield_lvmod.R


## 参考文献

ヌルクラインの説明については [Hirsch・Smale・Devaney 力学系入門（共立出版）](https://www.kyoritsu-pub.co.jp/book/b10003811.html) を参考にしました.

R のコードについては [StatModeling Memorandum - SIRモデルからはじめる微分方程式と離散時間確率過程（前編）](https://statmodeling.hatenablog.com/entry/sir-model-ode-1) を参考にしました.