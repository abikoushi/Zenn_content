---
title: "行列の指数関数と常微分方程式についてのイントロ"
emoji: "🪛"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [R, 微分方程式, 線形代数]
published: true
---
## はしがき

行列の指数関数と常微分方程式について入門するためのちょっとしたノートです．独学者が行列の指数関数にいきなり出会ってしまっても怯まないように，という願いが込められています．

先に主に参考にした文献を紹介しておきます．

1. [長谷川浩司『線形代数』](https://www.nippyo.co.jp/shop/book/6704.html)
2. [Hirsch・Smale・Devaney『力学系入門』](https://www.kyoritsu-pub.co.jp/book/b10003811.html)

タイトルの通り，1つめの方は線形代数が主眼で，2つめの方は力学系が主眼です．


## 行列の指数関数入門

### 線形1階常微分方程式

まず，自明でない中でおそらく一番かんたんな微分方程式である次を考える．

$$
\frac{d}{dt} X(t) = a X(t). \tag{1}
$$

$a$ は定数とした．この式は，

1. $a$ が正のときは「 $X(t)$ が $t$ に応じてどんどん大きくなる様子を表し，その $t$ 時点での変化量は $X(t)$ 自身に比例する」と見ることができる．
2. $a$ が負のときは「 $X(t)$ が $t$ に応じてだんだん小さくなる様子を表し，その $t$ 時点での変化量は $X(t)$ 自身に比例する」と見ることができる．


1 は，ネズミ算のように人口がどんどん増える様子のモデルに使えるかもしれない．2 はタンクから水圧で水が流出して水量がだんだん減っていく様子のモデルに使えるかもしれない．[『微分方程式で数学モデルを作ろう』](https://www.nippyo.co.jp/shop/book/1240.html) にはこのような例示がたくさん出ている．

このように微分方程式を「変化の仕方」に対して仮定をおくことに対応させると，それを解くことは「変化した結果どうなるか」に対応させられる．そのため，時間に応じて変化する現象の数理モデルでは，微分方程式は基本的な道具の一つとして広く使われている．

次に (1) の方程式を解くことを考える．ちなみに， 数学で「(1) の方程式を解く」というときは，特に断りがなければ「(1) の式が成り立つような $X(t)$ をすべて見つける」という意味であることが多い．

$X(t) = C \exp(a t)$ が (1) の式を満たすような関数であることは両辺を微分してみれば確かめられる．ここで $C$ は任意の定数（a.k.a. 積分定数）である．

実際，(1) は次のように書き換えることができるから，

$$
\frac{1}{X(t)}  dX(t)  = a dt 
$$

この両辺を積分して原始関数を求めると

$$
\log X(t)  = a t + C_0.
$$

ここで $C_0$ は積分定数である. 積分定数は任意なので $C=\exp(C_0)$ と改めておくと 

$$
X(t)  = C\exp(a t)
$$
 
が得られる．また $t=0$ の時点での値 $X(0)$ がなにか一つに決まっているとしたら $C=X(0)$ でなければならない．このような条件を **初期条件** と呼ぶ．

### 連立させてみる

いま，次の連立微分方程式を考えたい．

$$
\begin{aligned}
\frac{d}{dt} X(t) &= -a X(t) ,\\
\frac{d}{dt} Y(t) &= a X(t) - b Y(t).
\end{aligned} \tag{2}
$$

ひとまず， $a$ , $b$ は正の実数として， $X(t)$ が減った分 $Y(t)$ が増え， $Y(t)$ は $Y(t)$ 自身の量に比例して減る，という様子を表している．この微分方程式に実際に応用例があり得ることは，

- [Wikipedia - 反応速度式](https://ja.wikipedia.org/wiki/反応速度式) ；「連鎖反応」のところを参照
- [StatModeling Memorandum - 陽に解ける常微分方程式を使ったモデル](https://statmodeling.hatenablog.com/entry/ode-explicit)

などを見るとわかる．

(2) の組は次のように書くこともできる．

$$
\frac{d}{dt}
 \begin{pmatrix}X(t) \\
 Y(t)\end{pmatrix} = 
A \begin{pmatrix}X(t) \\
 Y(t)\end{pmatrix}. \tag{3}
$$

左辺の $d/dt$ は要素ごとの微分を表す．右辺の行列 $A$ は

$$
A = \begin{pmatrix} -a & 0 \\
 a & -b \end{pmatrix} 
$$

とおいた．(1) のときからの類推で (2) の解は

$$
 \begin{pmatrix}X(t) \\
 Y(t)\end{pmatrix} = 
\exp(At)
  \begin{pmatrix}X(0) \\
 Y(0)\end{pmatrix} \tag{4}
$$

とできそうな気がする．しかし，この時点では解 (4) はただの連想で論理的な考察ではないし，そもそも $\exp(At)$ がなんなのかよくわからないのでそれを考えていく．


突然で恐縮だが，行列 $A$ として固有値 $\lambda_i$ と固有ベクトル $\boldsymbol{v}_i$ で対角化できるものを考える．

$$
A = (\boldsymbol{v}_1, \boldsymbol{v}_2) \begin{pmatrix} \lambda_1 & 0 \\
0 & \lambda_2\end{pmatrix}   (\boldsymbol{v}_1, \boldsymbol{v}_2)^{-1}
$$

固有ベクトルを並べた行列を $P$, 固有値を対角成分に持つ行列を $\Lambda$ とおき，この式を $A=P \Lambda P^{-1}$ のようにも表すことにする．

さて $A^2=AA$ と書くことにすると，これについてはさほど違和感もないだろうし曖昧さもない． $A^2$ は次のようにも書ける．

$$
\begin{aligned}
A^2 & =P \Lambda P^{-1}\cdot P  \Lambda P^{-1}\\
&=P \Lambda  \Lambda P^{-1}\\
&=P \begin{pmatrix} \lambda_1^2 & 0\\ 
0 & \lambda_2^2\end{pmatrix}  P^{-1}.
\end{aligned}
$$

同様，2より大きい整数 $k$ に対しても，

$$
A^k  =P \begin{pmatrix} \lambda_1^k & 0 \\
0 & \lambda_2^k\end{pmatrix}  P^{-1}.
$$

これよりためしに行列の指数関数を

$$
\exp(At)  =P \begin{pmatrix} \exp(\lambda_1 t)  & 0 \\
0 & \exp(\lambda_2t)\end{pmatrix}  P^{-1}
$$

としてみる．これは成分ごとの微分に対し，

$$
\begin{aligned}
\frac{d}{dt} \exp(A t)  &= \frac{d}{dt}P \begin{pmatrix} \exp(\lambda_1 t) & 0 \\ 
  0 & \exp(\lambda_2 t)\end{pmatrix}  P^{-1} \\
&=P \begin{pmatrix} \lambda_1 \exp(\lambda_1 t) & 0 \\
  0 & \lambda_2 \exp(\lambda_2 t)\end{pmatrix}  P^{-1} \\
&= P \Lambda \begin{pmatrix} \ \exp(\lambda_1 t) & 0 \\
  0 &  \exp( t)\end{pmatrix}  P^{-1}\\
&=P \Lambda P^{-1} \cdot P\begin{pmatrix} \ \exp(\lambda_1 t) & 0 \\
  0 &  \exp(t)\end{pmatrix}  P^{-1}\\
&= A \exp(At)
\end{aligned}
$$

となる．これは望ましい性質である．なぜならこの関係を使って (4) を両辺微分してみると，(3)の式を満たしていることがわかるからである．

### 微分方程式の解の一意性

実は，これが(3)のすべての解である．(3)の任意の解を

$$
U(t) = \begin{pmatrix} u_1(t) \\ u_2(t) \end{pmatrix}
$$

とおく．$U(t)$ は (3) を満たすので

$$
\frac{d}{dt} U(t) = A U(t)
$$

である．

$$
\begin{aligned}
\frac{d}{dt} U(t) \exp(-At) &= A U(t) \exp(-At)-A U(t) \exp(-At)\\
&=0
\end{aligned}
$$

より，$U(t) \exp(-At)$ は定数ベクトルである．この定数ベクトルを $\boldsymbol{k}$ とおくと

$$
U(t) = \exp(At)\boldsymbol{k}
$$

と任意の解が (4) の形で表せることがわかる．

このように行列をつかって微分方程式が扱えそうなことがわかったことで，より多くの連立方程式をまとめて扱うなど，いろいろな可能性が広がる．

## R による数値的確認

微分方程式を数値的に解くには `deSolve` というパッケージが使える．(3)のような微分方程式を次のように書ける．

```r
library(deSolve)
modLinear <- function(Time, State, A) {
  return(list(A%*%State))
}
```


適当なパラメータを選び，解を求めてみる．

```r
a = 0.2
b = 0.1
A = matrix(c(-a, 0,
              a, -b), 2, 2, byrow = TRUE)
yini  <- c(X = 1, Y = 0)
times <- seq(0, 30, by = 0.5)
out   <- ode(yini, times, modLinear, A)
```

一方，(4)は下記の関数 `Solm` のように書ける．

```r
Solm = function(t, P, Pinv, values, yini){
  (P %*% diag(exp(t*values)) %*% Pinv) %*% yini  
}

ei_A = eigen(A)
P_A = ei_A$vectors
Pinv_A = solve(ei_A$vectors)

ts = 0:30
sol1 = sapply(ts, Solm, P=P_A, Pinv=Pinv_A, values=ei_A$values, yini=yini)
matplot(out_A[,1], out_A[,-1], type = "l", xlab = "t", ylab = "U(t)")
matpoints(ts, t(sol1))
```

![](/images/intro_matexp/sol1.png)
*数値解と解析解の比較*

数値解と解析解が一致していることがわかる．

## 固有値が虚数のとき

行列の固有値が実数の範囲に収まらないときも確認してみたい．

引き続き行列を使って次のように表せる微分方程式を扱う．

$$
\frac{d}{dt}U(t) = B U(t)  \tag{5}
$$

行列 $B$ として次のものを考えると，この固有値は $\pm i$ である（ここでの $i$ は虚数単位）．

$$
B=\begin{pmatrix}
0 & 1 \\
-1 & 0 \\
\end{pmatrix}.
$$

数値的に確認してみよう．

```r
B = matrix(c(0, 1,
             -1, 0), 2, 2, byrow = TRUE)
ei_B = eigen(B)
print(ei_B)
```

```
> ei_B
eigen() decomposition
$values
[1] 0+1i 0-1i

$vectors
                     [,1]                 [,2]
[1,] 0.7071068+0.0000000i 0.7071068+0.0000000i
[2,] 0.0000000+0.7071068i 0.0000000-0.7071068i

```

ちなみに `0.7071068` は $1/\sqrt{2}$ である．

先ほどと同様に数値解と行列の指数関数で表せる解析解を比べてみる．

```r
P_B = ei_B$vectors
Pinv_B = solve(P_B)
out_B   <- ode(yini, times, modLinear, B)

ts = 0:30
sol2 = sapply(ts, Solm, P=P_B, Pinv = Pinv_B, values=ei_B$values, yini=yini)
matplot(out_B[,1], out_B[,-1],type = "l", xlab = "t", ylab = "U(t)")
matpoints(ts, t(Re(sol2)))
```

![](/images/intro_matexp/sol2.png)
*数値解と解析解の比較*

数値解と解析解が一致していることがわかる．ただしここでは，虚部はほぼ0に近い値になっているので `Re()` を使って実部のみをプロットした．

```
> print(head(t(sol2)))
                         [,1]                     [,2]
[1,]  1.0000000+0.000000e+00i  0.0000000+1.014654e-17i
[2,]  0.5403023-1.377255e-17i -0.8414710-9.760998e-18i
[3,] -0.4161468-9.323820e-19i -0.9092974+2.995571e-18i
[4,] -0.9899925-2.393861e-18i -0.1411200+2.743162e-17i
[5,] -0.6536436-2.718334e-17i  0.7568025-2.447937e-17i
[6,]  0.2836622+1.062310e-17i  0.9589243+1.095779e-17i
```

図から，解は三角関数のように周期的な関数である．実は「三角関数のよう」というか三角関数そのものなのだが，ここではそれは置いておいて，この微分方程式がどういう状況を表すのに使えそうかだけ述べておく．

(5) の式はその要素を書き下すと次のようになる．

$$
\frac{d}{dt} 
\begin{pmatrix}
X(t)\\
Y(t)
\end{pmatrix}
=\begin{pmatrix}
Y(t)\\
-X(t)
\end{pmatrix}
$$

この式は， $Y$ が大きいときは $X$ の増加量も大きく，$Y$ の増加量は $X$ が少ないほど大きいことを表している．$X$ と$Y$ が捕食者と被食者のような関係にあるときの個体数の増減を表すのに使えるかもしれない．

また，冒頭で触れた [長谷川浩司『線形代数』](https://www.nippyo.co.jp/shop/book/6704.html) や [Hirsch・Smale・Devaney『力学系入門』](https://www.kyoritsu-pub.co.jp/book/b10003811.html)　では(5)と同様の方程式をバネの振動を例に解説している．高校などで物理を既習の方はこちらのほうがイメージしやすいかもしれない．

ちなみに，数値計算をいろいろやるときは今回やったように線形代数の教科書に出ているような式をそのまま打ち込むと遅いというのが結構あるあるで，例えば[伊理正夫・藤野和建『数値計算の常識](https://www.kyoritsu-pub.co.jp/book/b10011417.html) の第5章は「逆行列よさようなら」である．行列の指数関数については [expm](https://cran.r-project.org/web/packages/expm/index.html) というパッケージがあるのでよりより実践的にはこういうパッケージを使うのがいいと思う．

おしまい．