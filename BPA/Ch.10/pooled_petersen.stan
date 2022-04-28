data{
  int m;
  int c;
  int r;
}
parameters{
  real<lower=(c - r + m)> N;
}
model{
  //priors
  N ~ gamma(1E-06,1E-06);
  //Likelihood
  r ~ binomial(m, c/N);
}
