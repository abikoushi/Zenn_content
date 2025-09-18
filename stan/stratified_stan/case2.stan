//case.stan
data{
  int N;
  array[N] int<lower=0,upper=1> Y;
  array[N] int<lower=0,upper=1> Z;
  array[N] int<lower=0,upper=1> X;
  real<lower=0> alpha;
}
parameters{
    matrix<lower=0, upper=1>[2,2] Xi;
    vector<lower=0, upper=1>[2] phi;
    real<lower=0, upper=1> delta;
}
model{
  for(i in 1:N){
    Z[i] ~ bernoulli(phi[X[i]+1]);
    X[i] ~ bernoulli(delta);
    Y[i] ~ bernoulli(Xi[Z[i]+1, X[i]+1]);
  }
  to_vector(Xi) ~ beta(alpha, alpha);
  delta ~ beta(alpha, alpha);
  phi ~ beta(alpha, alpha);
}
generated quantities{
  int D1;
  int D2;
  real<lower=0, upper=1> gamma;
  vector<lower=0, upper=1>[2] psi;
  gamma = phi[1]*(1-delta) + phi[2]*delta;
  psi[1] = (1-phi[2])*delta/(1-gamma);
  psi[2] = phi[2]*delta/gamma;
  {
    int Zast = bernoulli_rng(gamma)+1;
    int Yast0 = bernoulli_rng(Xi[Zast,1]);
    int Yast1 = bernoulli_rng(Xi[Zast,2]);
    D1 = Yast1 - Yast0;
  }
  {
    int Yast0;
    int Yast1;
    int Zast0;
    int Zast1;
    Zast0 = bernoulli_rng(phi[1])+1;
    Zast1 = bernoulli_rng(phi[2])+1;
    Yast0 = bernoulli_rng(Xi[Zast0,1]);
    Yast1 = bernoulli_rng(Xi[Zast1,2]);
    D2 = Yast1 - Yast0;
  }
}
