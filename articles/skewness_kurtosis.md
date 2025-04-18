---
title: "Julia で学ぶ歪度・尖度"
emoji: "👻"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [Julia, 確率分布]
published: true
---

## あらまし

モーメント母関数の対数を取ったものをキュムラント母関数という．
キュムラント母関数のテイラー展開を用いると，中心極限定理による正規分布への収束の速さが歪度・尖度でほぼ決まることがわかる．

また，Julia の [Distibutions.jl](https://juliastats.org/Distributions.jl/stable/) パッケージでは確率分布の歪度・尖度の計算も実装されているので，これを使って数値的にもたしかめてみる．

## 手計算してみる

このパートは主に [竹内啓『数理統計学: データ解析の方法』（東洋経済新報社）](https://books.google.co.jp/books/about/数理統計学.html?id=fSo8DwAAQBAJ&redir_esc=y) を参考にしている（ただしまったく同じではない）．

### 準備

まず用語を確認する．

確率変数 $X$ の **モーメント母関数** $\phi_X(t)$ を

$$
\phi_X(t) = E[e^{tX}]
$$

と定義すると，テイラー展開

$$
\phi_X(t) = \phi_X  (0) + \phi_X' (0) t + \frac{\phi_X '' (0)}{2} t^2 + \frac{\phi_X^{(3)} (0)}{3!} t^3 + \frac{\phi_X^{(4)} (0)}{4!}t^4 + \cdots
$$

の $j$ 次の項の係数  $\phi_X^{(j)} (0)$ が $j$ 次のモーメント $\mu_j = E[X^j]$ になる．ただし $\phi_X  (0)  =  1$ である．

モーメント母関数の対数を取った $\psi_X (t) = \log \phi_X(t)$ を **キュムラント母関数** と呼ぶ．

キュムラント母関数のテイラー展開

$$
\psi_X(t) = \kappa_0 + \kappa_1 t + \frac{\kappa_2 }{2} t^2 + \frac{\kappa_3}{3!} t^3 + \frac{\kappa_4}{4!}t^4 + \cdots
$$

の $j$ 次の項の係数  $\kappa_j$ を $j$ 次の **キュムラント** と呼ぶ．

キュムラントとモーメントの間には次のような関係がある．

$$
\begin{aligned}
\kappa_0 &= \log \phi(0) = \log 1 = 0,\\
\kappa_1 &= \left. \frac{d}{dt} \log \phi(t) \right| _{t=0} = \phi'(0)/\phi(0) = \mu_1\\
\kappa_2 &= \frac{ \phi''(0)\phi'(0) - \{\phi'(0)\}^2}{\{\phi(0)\}^2} = \mu_2-\mu_1^2\\
&\vdots
\end{aligned}
$$

$\mu_2-\mu_1^2$ は $E[X^2] - \{E[X]\}^2$ なので $X$ の分布の分散である．$\sigma^2 = \mu_2-\mu_1^2$ と書くことにする．


### 中心極限定理と歪度・尖度

ここでは $j$ 次のモーメントがすべて存在する（積分が収束する）分布の場合のみを考える．

独立に同じ分布に従う確率変数 $X_i$ ($i=1,\ldots,n$) の標本平均 $\displaystyle \bar{X} = \frac{1}{n} \sum_{i=1}^n X_i$ が正規分布に収束することを示したい．

ところで，標準正規分布 $\mathcal{N}(0, 1)$ に従う確率変数 $Z$ を $Y=aZ+b$ （$a,b$は定数）と変数変換した $Y$ は平均 $\mu$，分散 $a^2$ の正規分布 $\mathcal{N}(b, a^2)$ に従う．

そのため，$\bar X$ が正規分布に従うことを示すために, 標本平均を $X_i$ の分布の平均 $\mu$ と標準偏差 $\sigma$ で変数変換 $Y= \sqrt{n}(\bar{X} - \mu)/\sigma$ が標準正規分布に従うことを示す．

$Y$ のモーメント母関数は，

$$
\begin{aligned}
\phi_Y(t) &= E[e^{t Y}] = E[\exp( \sqrt{n}(\bar{X} - \mu)/\sigma )]\\
&= E[\exp( \sqrt{n} \bar{X}/\sigma)]\exp( -\sqrt{n}\mu/\sigma) \\
&= \left(E\left[\exp \left( \frac{X_1}{\sigma\sqrt{n}}\right)\right] \right)^n\exp( -\sqrt{n}\mu/\sigma ).
\end{aligned}
$$

2行目から3行目の変形で独立同分布の仮定を使った．

$Y$ のキュムラント母関数をテイラー展開すると，「準備」で述べたキュムラントとモーメントの関係から，

$$
\begin{aligned}
\psi_Y(t) &= \log \phi_Y(t) \\
&= n\left\{\kappa_0 + \frac{\kappa_1}{\sqrt{n}\sigma} t + \frac{\kappa_2 /(\sqrt{n}\sigma)^2}{2} t^2 + \frac{\kappa_3/(\sqrt{n}\sigma)^3}{3!} t^3 + \frac{\kappa_4/(n^2 \sigma^4)}{4!}t^4 + \cdots \right\} -\sqrt{n}\mu/\sigma \\
&= \frac{\kappa_2 /\sigma^2}{2} t^2 + \frac{\kappa_3/\sigma^3}{\sqrt{n} \, 3!} t^3 + \frac{\kappa_4/\sigma^4}{n \, 4!}t^4 + \cdots
\end{aligned}
$$

と整理できる．右辺の第一項は標準正規分布のキュムラント母関数である．
キュムラント母関数は一意であり，他の項は $O(1/\sqrt{n})$ で 0に近づくので $n$ が大きいとき $Y$ が標準正規分布に近づくことがわかった．

$1/\sqrt{n}$ と同じくらいの速さで 0 に近づく3次の項の係数 $\kappa_3/\sigma^3$ を **歪度**（skewness），$1/n$ と同じくらいの速さで 0 に近づく4次の項の係数 $\kappa_4/\sigma^4$ を **尖度**（kurtosis）と呼ぶ．

歪度が 0 のときは正規分布への収束の速さが $O(1/n)$ になって，その速さが尖度でほぼ決まるということである．

中心極限定理は数理統計学の教科書にはだいたい出ていると思うが，証明はモーメント母関数を使うもの（例えば，野田一夫・宮岡悦良『入門演習　数理統計』共立出版）や，特性関数を使うもの（例えば，久保川達也『現代数理統計学の基礎』共立出版）が多い印象がある．

ここで紹介した [竹内啓『数理統計学: データ解析の方法』（東洋経済新報社）](https://books.google.co.jp/books/about/数理統計学.html?id=fSo8DwAAQBAJ&redir_esc=y) のようにキュムラント母関数を使う方針だと，計算の難しさはほぼ同じのまま $n \to \infty$ での収束先だけでなく収束の速さも知れるという良さがある．このメリットはけっこう大きいと思う．

ただし特性関数を使うとモーメント母関数が存在しない場合についても示せるというメリットがある．

## 数値計算してみる

また，Julia の [Distibutions.jl](https://juliastats.org/Distributions.jl/stable/) パッケージでは確率分布の歪度は `skewness` という関数で，尖度は `kurtosis` という関数で実装されている．

こんなふうにすると正規分布の歪度・尖度が出る：

```julia
using Distributions
using Random
using Statistics
using Plots

N0 = Normal(0,1)
skewness(N0) #0
kurtosis(N0) #0
```

歪度が 0 の分布には例えば一様分布やロジスティック分布がある．こんな形をしている．

```julia
U1 = Uniform(-1,1)
skewness(U1) #0.0
kurtosis(U1) #-1.2

L1 = Logistic(0,1)
skewness(L1) #0.0
kurtosis(L1) #1.2

p = plot(x -> pdf(N0,x), -4, 4, tick_direction=:out, label="Normal", linestyle=:solid)
plot!(p, x -> pdf(U1,x), label="Uniform", linestyle=:dash)
plot!(p, x -> pdf(L1,x), label="Logistic", linestyle=:dot)
```

![](/images/skewness_kurtosis/density1.png)

次のように標本分布 $Y= \sqrt{n}(\bar{X} - \mu)/\sigma$ を繰り返し求める関数を宣言する．

```julia
function simCLT(dist, size, iter)
    out = zeros(iter)
    rng = Random.default_rng()
    for i in 1:iter
        X = rand(rng, dist, size)
        out[i] = (mean(X) - mean(dist)) * sqrt(size)/std(dist)
    end
     return out
end
```

一様分布のときの $Y$ のヒストグラムを見ると $n=10$ でも結構正規分布に近い．青い線で標準正規分布の密度関数を重ねてある．

```julia
col1 = palette(:default)[1:3]
outU1 = simCLT(U1, 10, 10000)
h = histogram(outU1, legend=false, normalize=:pdf, color="lightgray", tick_direction=:out)
plot!(h, x->pdf(Normal(),x), linewidth=2, color=col1[1])
```

![](/images/skewness_kurtosis/hist1.png)


ロジスティック分布のときも同様である．

```r
outL1 = simCLT(L1, 10, 10000)
h = histogram(outL1, legend=false, normalize=:pdf, color="lightgray", tick_direction=:out)
plot!(h, x->pdf(Normal(),x),  linewidth=2, color=col1[1])
```

![](/images/skewness_kurtosis/hist2.png)

一方，歪度が 0 より大きい分布の一例として指数分布をみる．

```julia
E1 = Exponential(1)

skewness(E1) #2.0
kurtosis(E1) #6.0

p = plot(x -> pdf(N0,x), -4, 4, tick_direction=:out, label="Normal", linestyle=:solid)
plot!(x->pdf(E1,x), 0,4, label="Exponential", linestyle=:dash,tick_direction=:out)
```

![](/images/skewness_kurtosis/density2.png)


指数分布のときの $Y$ のヒストグラムを見ると $n=10$ では正規分布から結構離れている．青い線で標準正規分布の密度関数を重ねてある．

```julia
outE11 = simCLT(E1,  10, 10000)
h = histogram(outE11, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color=col1[1])
```

![](/images/skewness_kurtosis/hist3.png)

$n=100 ~(=10^2)$ くらいまで大きくするとちょっと近づいてくる．

```julia
outE12 = simCLT(E1, 100, 10000)
h = histogram(outE12, legend=false, normalize=:pdf, color="lightgray")
plot!(h, x->pdf(Normal(),x), color=col1[1])
```

![](/images/skewness_kurtosis/hist4.png)

（おしまい）