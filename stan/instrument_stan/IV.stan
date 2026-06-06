data {
  int<lower=0> N;
  vector[N] Y;
  vector[N] X;
  vector[N] Z;
}
parameters {
  vector[N] U;
  real alpha_xz;
  real alpha_yx;
  real<lower=0> alpha_xu;
  real<lower=0> alpha_yu;
}
model {
  for(i in 1:N){
    Z[i] ~ normal(0,1);
    U[i] ~ normal(0,1);
    X[i] ~ normal(alpha_xz*Z[i] + alpha_xu*U[i], 1);
    Y[i] ~ normal(alpha_yx*X[i] + alpha_yu*U[i], 1);
  }
  alpha_xz ~ normal(0, 10);
  alpha_xu ~ normal(0, 10);
  alpha_yx ~ normal(0, 10);
  alpha_yu ~ normal(0, 10);
}
generated quantities {
  real TE;
  {
  real Y1;
  real Y0;
  real Uast = normal_rng(0, 1);
  Y1 = normal_rng(alpha_yx+alpha_yu*Uast, 1);
  Y0 = normal_rng(alpha_yu*Uast, 1);
  TE = Y1-Y0;
  }
}
