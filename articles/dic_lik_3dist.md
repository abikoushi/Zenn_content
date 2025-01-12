---
title: "2重区間打ち切りされた観測についての尤度"
emoji: "🔖"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 統計]
published: true
---

## あらまし

[Reich et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19598148/) は感染症の潜伏期間を推定する際によく参照される文献で2重区間打ち切り（doubly interval censored）のある観測について論じている．

Reich さんたちは Reich et al. (2009) の考察に基づく R のパッケージ `coaseDataTools` も提供しており，このパッケージでは2重区間打ち切りされた観測に対してワイブル分布，ガンマ分布，対数正規分布などを当てはめることができる．

`coaseDataTools` の中身を見てみると尤度の中の積分は数値積分で評価している様子だったが，[ウルフラム・アルファ](https://ja.wolframalpha.com) を使って計算してみたらワイブル分布，ガンマ分布，対数正規分布についてはこの積分が閉じた形（といっても特殊関数は使う）で書けたのでここで報告する．

ガンマ関数などの特殊関数はそれ用の数値計算のプログラムがうまく作られていたりするので特殊関数を使って書けると計算が速くなることある（もともとそこまで時間のかかる計算ではなかったので計算を速くしたいという意欲はそんなにはないかもしれないが……）．


## 2重区間打ち切り

準備としてまず観測の枠組みを説明する．

感染の起こった時点に対応する確率変数を $E$, 発症した時点に対応する確率変数を $S$ とする.

潜伏期間は $Y=S-E$ であり，$Y$ の分布 $F(y)$ が知りたい．

このとき尤度は $F(y)$ の密度関数 $f(y)$ と $E$, $S$ の実現値 $e$, $s$ を用いて,

$$
\ell = f(s-e) = f(s-e)
$$

と書ける．

ここで $S$ の正確な時点はわからず $[S_L, S_R]$ として観測されたとしよう．

このようにイベントが起こったことはわかっているがイベントの時点が幅を持って観測される状況を区間打ち切りと呼ぶ．

この場合, $S$ が $[S_L, S_R]$ という期間の中で定常的であるとすると（あるいは $S$ に対して flat prior $\propto 1$ を用いることにすると）尤度は

$$
\begin{aligned}
\ell &= \int_{S_L}^{S_R} f(s-e) \,de \\
&= F(S_R-e)-F(S_L-e)
\end{aligned}
$$

と書ける．

さらに，感染した時点 $E$ が特定されることはまれであり，$E$ は感染の流行していた場所に滞在していた期間など，なんらかのリスク因子に曝露していた期間 $[E_L, E_R]$ として観測されることが多いと考えられる．

[Reich et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19598148/) はこのような観測を2重区間打ち切り（doubly interval censored）と呼んだ．

$E$ もまた $[E_L, E_R]$ という期間の中で定常的であるとすると尤度は，

$$
\begin{aligned}
\ell &=\int_{E_L}^{E_R}  \int_{S_L}^{S_R}f(s-e) \,de \, ds\\
&= \int_{E_L}^{E_R} F(S_R-e)-F(S_L-e) \, de
\end{aligned}
$$

と書ける. この記事ではこの尤度を具体的に計算してみる．

（ところで，感染の打ち切りは発生するだろうなと思うけど発症した日は大体わかりそうな気がする．現実のどういう場面で発症日が打ち切りになるのかいまいちわかってないのでだれか教えてください．）

## 準備：$1-F(y)$ の積分

先の積分について具体的に計算する前に，少しだけ記号を準備しておく．

分布 $F(y)$ の平均を $\mu$ として（平均の存在は仮定），$F(y)$ が 0 でない値をとる範囲（これを関数の台とかサポートとかいう）を $[0,\infty)$ とする（過去にさかのぼって発症するようなことは考えない）．

あとちょっと面倒を避けるために $F(y)$ としてとりあえず連続型の分布だけを考えることにする．

打ち切りデータを扱う生存時間分析の文脈だと生存関数 $S(y)=1-F(y)$ を考えることが多いので，ここでもそれにならって $1-F(y)$ の積分をまず考えることにする．

$$
F(S_R-e)-F(S_L-e) = 1-F(S_R-e)-(1-F(S_L-e))
$$

であり，$1-F(y)$ の積分がわかると $F(y)$ の積分もわかる．そして，いま

$$
\mu = \int_0^{\infty} (1-F(y)) \, dy
$$

が成り立つ．なぜかというと，部分積分により, 

$$
\begin{aligned}
\mu &= \int^{\infty}_{0}y f(y) \, dy  \\
&= \int^{\infty}_{0}y  \{-(1 - F(x)) \}' \, dy \\
&= \biggl[ y  \{-(1 - F(y)) \}\biggr]^{\infty}_{0} -  \int^{\infty}_{0} \{-[1 - F(x)] \} \, dx \\
&= \int^{\infty}_{0} [1 - F(x) ]\, dx 
\end{aligned}
$$

が成り立つからである．最後の変形で第1項が 0 になるのは，

$$
y  \{1 - F(y) \}= y \int_{y}^{\infty} f(u) \, du \le \int_{y}^{\infty}u  f(u) \, du 
$$

であり，

$$
\lim_{u \to \infty} \int_{u}^{\infty}y  f(y) \, dy = 0
$$

であることからわかる. 

ここから生存関数を $\mu$ で割って正規化し，

$$
f^{eq}(y) = \frac{1-F(y)}{\mu}
$$

とすると，$f^{eq}(y)$ も全区間で積分して1となるから確率密度関数になっている．

この確率密度を持つ分布を [Ross (1995) "Stochastic Processes"](https://www.wiley.com/en-us/Stochastic+Processes%2C+2nd+Edition-p-9780471120629) にならい，均衡分布（equibillium distribution）と呼ぶことにする．

## 均衡分布

ガンマ分布，ワイブル分布，対数正規分布について，均衡分布の分布関数を求めてみたところ以下に示すような形になった．

興味のある方は[ウルフラム・アルファ](https://ja.wolframalpha.com) を開いて `integral[ Gamma(a, x/b) , x]` などと打ち込んでみてほしい．

分布関数を使って書けたので（この記事では扱わないが）[Stan](https://mc-stan.org) などの確率的プログラミング言語でも実装しやすいのではないかと思う．

### ガンマ分布の場合

$\bar \gamma(a,t)$, $\bar \Gamma(a,t)$ を全区間での積分が1になるよう正規化された不完全ガンマ関数，

$$
\bar \gamma (a,x)= \frac{1}{\Gamma(a)} \int _{0}^{x}t^{a-1}\,e^{-t}\,dt
$$

$$
\bar \Gamma (a,x)= \frac{1}{\Gamma(a)} \int _{x}^{\infty}t^{a-1}\,e^{-t}\,dt
$$

とする. 形状パラメータ $a$，尺度パラメータ $b$ のガンマ分布の均衡分布の分布関数は，

$$
F^{eq}(y) = \bar \gamma(a+1, x/b)+ \frac{x}{ab} \bar \Gamma(a, x/b).
$$

である．

正規化されたガンマ関数はガンマ分布の分布関数で尺度パラメータ $b=1$ のときと同じであるから，R で実装するときは次のように書ける：

```R
eqgamma <- function(x, shape, scale){
  pgamma(x/scale, shape+1)+(x/scale)*pgamma(x/scale, shape, lower.tail=FALSE)/shape
}
```

### ワイブル分布の場合

形状パラメータ $a$，尺度パラメータ $b$ のワイブル分布の均衡分布の分布関数は，

$$
F^{eq} (y) = \bar \gamma(1/a, (x/b)^a).
$$

R での実装例は：

```R
eqweibull <- function(x, shape, scale){
  pgamma((x/scale)^(shape), 1/shape)
}
```

### 対数正規分布

$\operatorname{erf}(x)$, $\operatorname{erfc}(x)$ をそれぞれ誤差関数と補誤差関数，

$$
\operatorname {erf} \left(x\right)={\frac {2}{\sqrt {\pi }}}\int _{0}^{x}e^{-t^{2}}\,dt,
$$

$$
\operatorname {erfc} \left(x\right)={\frac {2}{\sqrt {\pi }}}\int _{x}^{\infty}e^{-t^{2}}\,dt
$$

とすると，対数平均パラメータ $\mu$, 対数標準偏差パラメータ $\sigma$ の対数正規分布の均衡分布の分布関数は，

$$
\begin{aligned}
F^{eq} (y) = &\frac{1}{2} \left(1+x \operatorname{erfc}\left(\frac{\log(x) - \mu}{\sigma \sqrt(2)} \right)\exp(-(\mu + (\sigma^2)/2)\right) \\
&- \frac{\operatorname{erf}(\mu + \sigma^2 - \log(x)}{\sigma \sqrt(2))}.
\end{aligned}
$$

R での実装例は：

```r
erf <- function(x){
  2*pnorm(sqrt(2)*x) - 1
}

erfc <- function(x){
  2*pnorm(-sqrt(2)*x) 
}
eqlnorm <- function(x, meanlog, sdlog){
  ifelse(
    x < 0, 
    0,
    0.5*(1+x*erfc((log(x) - meanlog)/(sdlog*sqrt(2)))*exp(-(meanlog + (sdlog^2)/2)) - erf((meanlog + sdlog^2 - log(x))/(sdlog * sqrt(2))))
  )
}
```


## 2重区間打ち切りされた観測についての尤度

再び2重区間打ち切りされた観測についての尤度に戻る．

これまでの結果を使って，2重区間打ち切りされた観測についての尤度を均衡分布で表してみる．

積分区間に注意して3つの場合にわける．

![場合分けの図](/images/dic_lik_3dist/threecases.png)


$E_R \le S_L$ のとき,

$$
\ell = \mu (F^{eq}(S_R-E_R)-F^{eq}(S_L-E_R)-F^{eq}(S_R-E_L)+F^{eq}(S_L-E_L))
$$

$LS < RE$ のとき,

$$
\ell = \mu F^{eq}(S_R-E_R) + E_R + S_L - \mu \{ F^{eq}(S_R-E_L)+F^{eq}(S_L-E_L) \}.
$$


$RS = RE$ のとき,

$$
\ell = (RS-LS) - \mu(F^{eq}(RS-LE)-F^{eq}(LS-LE)) . 
$$

感染が起こった時点は少なくとも発症よりは前のはずなので，$RS \le RE$ の場合は $RS = RE$ と考える．

独立同分布の仮定のもとで $n$ 件の観測が得られたとき，最尤推定のために最大化すべき尤度はこれらの因子 $n$ 個をすべてかけ合わせたものである．


## `coaseDataTools` との比較

`coaseDataTools` パッケージのデモ用に入っているインフルエンザの潜伏期間のデータで，最尤推定の計算時間を比較してみる．

R のコードは少し長くなるので最後にまとめおく．

期待通り速くなり，結果はおおよそ一致した．

```r
#上が coaseDataTools で実装されている尤度, 下は僕の実装
#ガンマ分布
> print(bm_g)
                                                                    test
1 flu_opt_g0 <- optim(c(0, 1), loglikhd, dat = fluA.inc.per, dist = "G")
2                   flu_opt_g <- optim(c(0, 1), lpG, dat = fluA.inc.per)
  replications elapsed relative user.self sys.self user.child sys.child
1          100  36.063    5.408    35.949    0.117          0         0
2          100   6.668    1.000     6.510    0.027          0         0
> print(all.equal(flu_opt_g0$par, flu_opt_g$par))
[1] TRUE

#ワイブル分布
> print(bm_w)
                                                                    test
1 flu_opt_w0 <- optim(c(0, 1), loglikhd, dat = fluA.inc.per, dist = "W")
2                   flu_opt_w <- optim(c(0, 1), lpW, dat = fluA.inc.per)
  replications elapsed relative user.self sys.self user.child sys.child
1          100  31.457    7.398    31.174    0.160          0         0
2          100   4.252    1.000     4.238    0.011          0         0
> print(all.equal(flu_opt_w0$par, flu_opt_w$par))
[1] TRUE

#対数正規分布
> print(bm_l)
                                                                    test
1 flu_opt_l0 <- optim(c(0, 1), loglikhd, dat = fluA.inc.per, dist = "L")
2                   flu_opt_l <- optim(c(0, 1), lpL, dat = fluA.inc.per)
  replications elapsed relative user.self sys.self user.child sys.child
1          100  52.890    5.015    52.598    0.206          0         0
2          100  10.547    1.000    10.526    0.024          0         0
> print(all.equal(flu_opt_l0$par, flu_opt_l$par))
[1] TRUE
```

せっかくなので生存関数もプロットしておく．

```r
c3 <- c("royalblue", "orangered", "forestgreen")
curve(pweibull(x, exp(flu_opt_w$par[1]), exp(flu_opt_w$par[2]), lower.tail = FALSE),
      xlim=c(0,5), col=c3[1], lty=2, ylab="CCDF", xlab="incubation period")
curve(pgamma(x, exp(flu_opt_g$par[1]), scale = exp(flu_opt_g$par[2]), lower.tail = FALSE),
      add=TRUE, col=c3[2], lty=3)
curve(plnorm(x, flu_opt_l$par[1], sdlog = exp(flu_opt_l$par[2]), lower.tail = FALSE),
      add=TRUE, col=c3[3], lty=4)
legend("topright",  legend = c("Weibull", "gamma", "lognormal"), lty=2:4, col=c3, lwd=1)

```

![生存関数3つ](/images/dic_lik_3dist/fluA.png)

選んだモデルの範囲だとどれもそれほど大きな違いはない．

横軸の単位は日だと思うが，ウイルスは変異とかもあるしこの潜伏期間が文字通りに解釈できるかはわからない（念のための注意）．

以下，R のコード：

```r
#install.packages("coarseDataTools")
library(coarseDataTools)
library(dplyr)
library(ggplot2)
library(rbenchmark)
#gamma dist
eqgamma <- function(x, shape, scale){
  pgamma(x/scale, shape+1)+(x/scale)*pgamma(x/scale, shape, lower.tail=FALSE)/shape
}
meangamma <- function(shape, scale){
  scale*shape
}

#weibull dist
eqweibull <- function(x, shape, scale){
  pgamma((x/scale)^(shape), 1/shape)
}
meanweibull <- function(shape, scale){
  scale*gamma(1+1/shape)
}

#lognormal dist
erf <- function(x){
  2*pnorm(sqrt(2)*x) - 1
}

erfc <- function(x){
  2*pnorm(-sqrt(2)*x) 
}
eqlnorm <- function(x, meanlog, sdlog){
  ifelse(
    x < 0, 
    0,
    0.5*(1+x*erfc((log(x) - meanlog)/(sdlog*sqrt(2)))*exp(-(meanlog + (sdlog^2)/2)) - erf((meanlog + sdlog^2 - log(x))/(sdlog * sqrt(2))))
  )
}
meanlnorm <- function(mu, sigma){
  exp(mu+sigma^2/2)
}


#doubly interval censored
evaluatelp0 <- function(EL, ER, SL, SR, par, eqcdf, meanf){
  ll = 0 
  mu = meanf(par)
  if(ER <= SL){
    ll = log(eqcdf(SR-ER, par) - eqcdf(SL-ER, par) - (eqcdf(SR-EL, par) - eqcdf(SL-EL, par))) + log(mu)
  }else if (SL < ER & ER < SR){
    ll = log((ER - SL) + mu*(eqcdf(SR-ER, par) - (eqcdf(SR-EL, par) - eqcdf(SL-EL, par))))
  }else if(SR <= ER){
    ll = log((SR - SL) - mu*(eqcdf(SR-EL, par) - eqcdf(SL-EL, par)) )  
  }
  return(ll)
}

# single interval
evaluatelp1 <- function(EL, ER, SL, SR, par, cdf){
  ll = log(cdf(SR-EL, par) - cdf(SL-EL, par))
  return(ll)
}

# no censored
evaluatelp2 <- function(EL, ER, SL, SR, par, logpdf){
  logpdf(SL-EL, par)
}

eval_lp <- function(par, EL, ER, SL, SR, type, dist){
  if(dist=="weibull"){
    eqcdf <- function(x, par){
      eqweibull(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    cdf <- function(x, par){
      pweibull(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    logpdf<- function(x, par){
      dweibull(x, shape = exp(par[1]), scale = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meanweibull(exp(par[1]), exp(par[2]))
    }
  }else if(dist == "gamma"){
    eqcdf <- function(x, par){
      eqgamma(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    cdf <- function(x, par){
      pgamma(x, shape = exp(par[1]), scale = exp(par[2]))
    }
    logpdf<- function(x, par){
      dgamma(x, shape = exp(par[1]), scale = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meangamma(exp(par[1]), exp(par[2]))
    }
  }else if(dist == "lognormal"){
    eqcdf <- function(x, par){
      eqlnorm(x, meanlog = par[1], sdlog = exp(par[2]))
    }
    cdf <- function(x, par){
      plnorm(x,  meanlog = par[1], sdlog = exp(par[2]))
    }
    logpdf<- function(x, par){
      dlnorm(x, meanlog = par[1], sdlog = exp(par[2]), log=TRUE)
    }
    meanf <- function(par){
      meanlnorm(par[1], exp(par[2]))
    }
  }else{
    warning("this distribution is not implemented")
  }
  N <- length(EL)
  ll = 0
  for(i in 1:N){
    if(type[i]==0){
      ll = ll + evaluatelp0(EL[i], ER[i], SL[i], SR[i], par, eqcdf = eqcdf, meanf = meanf)
    }else if(type[i]==1){
      ll = ll + evaluatelp1(EL[i], ER[i], SL[i], SR[i], par, cdf = cdf)
    }else if(type[i]==2){
      ll = ll + evaluatelp2(EL[i], ER[i], SL[i], SR[i], par, logpdf = logpdf)
    }else{
      warning("this censor type is not implemented")
    }    
  }
  return(ll)
}

data(fluA.inc.per)
#head(fluA.inc.per)


lpW <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "weibull"))  
}

lpG <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "gamma")) 
}

lpL <- function(par, dat){
  with(dat, -eval_lp(par, EL, ER, SL, SR, type = type, dist = "lognormal")) 
}

bm_w <- benchmark(flu_opt_w0 <- optim(c(0,1), loglikhd, dat=fluA.inc.per, dist = "W"),
  flu_opt_w <- optim(c(0,1), lpW, dat=fluA.inc.per), order = NULL)

print(bm_w)
print(all.equal(exp(flu_opt_w0$par),exp(flu_opt_w$par)))

bm_g <- benchmark(flu_opt_g0 <- optim(c(0,1),loglikhd,dat=fluA.inc.per,dist = "G"),
                  flu_opt_g <- optim(c(0,1), lpG, dat=fluA.inc.per), order = NULL)
print(bm_g)
print(all.equal(exp(flu_opt_g0$par),exp(flu_opt_g$par)))

bm_l <- benchmark(flu_opt_l0 <- optim(c(0,1),loglikhd, dat=fluA.inc.per,dist = "L"),
                  flu_opt_l <- optim(c(0,1), lpL, dat=fluA.inc.per), order = NULL)

print(bm_l)
print(all.equal(exp(flu_opt_g0$par),
                exp(flu_opt_g$par)))
```