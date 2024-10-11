---
title: "事後分布を正規分布で近似する変分ベイズ法"
emoji: "📑"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: []
published: false
---



[須山『ベイズ推論による機械学習入門』](https://amzn.to/3BHWbv1)ではロジスティック回帰とニューラルネットのところで近似事後分布として独立な正規分布を仮定して、変分パラメータを勾配法で推定するアルゴリズムが出てくる。

勾配を評価するときには再パラメータ化トリック（re-parameterization trick）というアイデアを使う。

これはニューラルネットワークに限らずどんな分布でも使えるアイデアで便利だと思ったので、まとめ直してみる。

** アルゴリズム

ここではパラメータの添字は省略する。

まずすべてのパラメータ $w$ に対して、平均0、分散 $\lambda^{-1}$ の正規事前分布 $p(w)$ を設定する。

[tex:\displaystyle \log p(w) =\frac{1}{2} \log \lambda - \lambda \frac{w^2}{2}+\mathrm{Const.}]

すべてのパラメータの近似事後分布として、平均[tex:\mu]、分散 [tex:\sigma^{2}] の正規分布 [tex:q(w)] を設定する。

[tex:\displaystyle \log q(w) = - \log \sigma -  \frac{(w-\mu)^2}{2\sigma^2}+\mathrm{Const.}]

標準偏差は正の数なので $\sigma = \exp(\rho)$ とおいて、$\rho$ を最適化する。

『ベイズ推論による機械学習入門』では、$\sigma = \log(1+\exp(\rho))$ となっているが、微分をかんたんにするため、指数関数にしてみた。

真の事後分布とのカルバック・ライブラ距離が近い正規分布を求めるのがアルゴリズムの目的。

標準正規乱数 [tex:\varepsilon] を使ってサンプル[tex: \tilde{w} = \mu + \sigma \tilde{\varepsilon}] を得るとこれは、近似事後分布からのサンプルそのものになる。観測された目的変数を [tex: Y]、 説明変数を[tex: X] として、尤度を [tex: p(Y|X,w)] と書くと真の事後分布と近似事後布のカルバックライブラ距離は、

[tex: \displaystyle \operatorname{KL}(q(w)\|p(w|Y,X))\\
 \displaystyle \approx \log q(\tilde w)- \log p(\tilde w) - \log p(Y|X,\tilde w)+\mathrm{Const.}]

と近似できる。

1個のサンプルの標本平均で期待値を近似しているわけで、モンテカルロEMとかを使ったことがある人は変に感じるかもしれないけれでも、これでもけっこううまくいく。もちろん複数サンプリングして平均をとってもよいが計算が大変になる。

あとはこれを微分して、勾配法で [tex:\mu] と [tex:\sigma] を更新してやる。

$$
 \log q(\tilde w) = - \log \sigma -  \frac{(\tilde w - \mu)^2}{2\sigma^2}+\mathrm{Const.}\\
= - \log \sigma -  \frac{(\mu+\sigma \varepsilon-\mu)^2}{2\sigma^2}+\mathrm{Const.}\\
 = - \log \sigma -  \frac{ \varepsilon^2}{2}+\mathrm{Const.}
$$

なので、
[tex: \displaystyle \frac{d}{d\mu}\log q(\tilde w) =0]
[tex: \displaystyle \frac{d}{d\sigma}\log q(\tilde w) = (-1/\sigma)]
[tex: \displaystyle \frac{d}{d\mu}\log p(\tilde w) =-\lambda \tilde w]
[tex: \displaystyle \frac{d}{d\sigma}\log p(\tilde w) = -\lambda \tilde w \varepsilon]

また、合成関数の微分なので、
[tex: \displaystyle \frac{d}{d\rho} \exp(\rho) =\exp(\rho)]
を忘れないようにかけてやる。 

まとめると次のアルゴリズムが得られる。
+ 学習率 [tex:\alpha] と 事前分布の精度パラメータ [tex: \lambda] を設定。
+ [tex:\mu] と [tex:\rho] を適当に初期化し以下を繰り返す。
+ すべてのパラメータに対して、標準正規乱数 [tex:\tilde{\varepsilon}] を使ってサンプル [tex: \tilde{w} = \mu + \exp(\rho) \tilde{\varepsilon}] を得る。 
+ [tex: g = -\frac{d}{dw} \log p(Y|X, \tilde w)] を計算。
+ [tex: g_\mu =  g + \lambda\tilde{w}] とする。
+ [tex: g_\rho = (g \tilde{\varepsilon} + \lambda \tilde{w} \tilde{\varepsilon})\exp(\rho) -1] とする。
+ [tex: \mu \leftarrow \mu + \alpha g_\mu] で更新
+ [tex: \rho \leftarrow \rho + \alpha g_\rho] で更新

** Julia のコード

Julia には自動微分のパッケージ ForwardDiff があるので対数尤度の微分のところは自分で計算しなくてよい（場合もある）。

今回は乱数で適当に作ったデータでポアソン回帰をやってみる。

こんなふうだ：

>|julia|

using Distributions
using ForwardDiff
using Random
using Plots
using Optim
using LinearAlgebra

function　poisonloss(beta,y,X)
    Xbeta = X*beta
    lambda = exp.(Xbeta)
    return -sum(y .* Xbeta - lambda)
end

function VIGD(f, par0, lambda, lr, maxiter::Int) 
    rng = Random.default_rng()
    len = length(par0)
    mu = randn(rng,len)
    rho = randn(rng,len)
    g(beta) = ForwardDiff.gradient(beta0 -> f(beta0),beta)
    logloss = zeros(maxiter)
    for i in 1:maxiter
    sigma = exp.(rho)
    epsilon = randn(rng,len)
    beta = mu + sigma.*epsilon
    fx = f(beta)
    gvec = g(beta)
    g_mu = gvec + lambda*beta
    g_rho = (gvec.*epsilon + lambda*beta.*epsilon).*sigma .- 1.0
    mu = mu - lr * g_mu
    rho = rho - lr * g_rho
    logloss[i] = fx
    end
    return mu, rho, logloss
end

rng = Random.default_rng()
Random.seed!(1)
x = sort(randn(rng,100))
X = [ones(100) x]
beta = [2.0,-1.0]

y = rand.(rng,Poisson.(exp.(X*beta)))

f(beta) = poisonloss(beta,y,X)
    
betaini = [0.0,0.0]
@time μ, ρ, logloss = VIGD(f, betaini, 0.0, 1.0e-4, 5000)
#計算時間は0.2秒くらい

plot(logloss,legend=false)

ε = randn(rng,2,1000)
betasmp = μ .+ exp.(ρ).*ε 

post = exp.(X*betasmp)
pred = rand.(rng,Poisson.(post))
predmean = mean(pred,dims=2)
lwr = [quantile(pred[i,:],0.025) for i in 1:100]
upr = [quantile(pred[i,:],0.975) for i in 1:100]

scatter(x,y,legend=false)
plot!(x,predmean, color="blue")
plot!(x,lwr,color="blue", linestyle = :dash)
plot!(x,upr,color="blue", linestyle = :dash)
png("./Desktop/plot.png")

opt = optimize(f,betaini,method=BFGS(),autodiff=:forward)
β = Optim.minimizer(opt)
se = sqrt.(diag(inv(Symmetric(ForwardDiff.hessian(f,β)))))
||<

poisonloss のところを書き換えればロジスティック回帰でもワイブル回帰でもできるはず。

予測分布の95%区間を点線でプロットした。

[f:id:abrahamcow:20201229143541p:plain]

計算が進むにつれて負の対数尤度が小さくなっていく様子です。

[f:id:abrahamcow:20201229143941p:plain]

ちなみに最尤推定の標準誤差にくらべると変分近似した事後分布の標準偏差はやや小さめに求まる。

>|julia|
julia> println(μ)
[1.999546299946068, -1.0194307666319837]

julia> println(β)
[2.0010916047667635, -1.0185383086714856]

julia> println(exp.(ρ))
[0.035814747235478735, 0.026642018884133416]

julia> println(se)
[0.04194236834667879, 0.031554017967900964]
||<

とはいえ、平均はともかく標準偏差のほうはイテレーションの回数によってけっこう変わる。

** もっと学びたい人のために

[https://arxiv.org/abs/1505.05424:title]

関連エントリ：[https://selfboast.hatenablog.jp/entry/2022/01/01/110919:title]
