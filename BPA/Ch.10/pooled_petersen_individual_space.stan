data {
  int<lower=0> M; // Size of augumented data set
  int<lower=0> T; // Number of sampling occasions
  int<lower=0> C; // Size of observed data set
  int<lower=0, upper=1> y[M, T]; // Capture-history matrix
  vector<lower=-6, upper=6>[C] bsize; // Body size
  int S; //number of reaches
  int<lower=0> reach[C]; // reach
}
transformed data {
  int<lower=0> s[M]; // Totals in each row

  for (i in 1 : M) {
    s[i] = sum(y[i]);
  }
}
parameters {
  real<lower=0, upper=1> omega; // Inclusion probability
  real<lower=0, upper=1> mean_p; // Mean detection probability
  real beta;
  vector[S] b2;
  real mu_size;
  real<lower=0> sd_size;
  vector<lower=-6, upper=6>[M - C] bsize_mis; // Missing data
  simplex[S] p_S;//proportion of fish tagged in each reach
}
transformed parameters {
}
model {
  matrix[M, T] logit_p;
  // Priors
  //  omega ~ uniform(0, 1);
  mean_p ~ uniform(0, 1);
  beta ~ std_normal();
  mu_size ~ std_normal();
  sd_size ~ std_normal();


  for (i in 1 : C) {
    logit_p[i,] = to_row_vector(rep_vector(logit(mean_p) + beta * bsize[i] + b2[reach[i]], C));
  }

  // Likelihood
  for (i in 1 : C) {
    bsize[i] ~ normal(mu_size, sd_size);
    reach[i] ~ categorical(p_S);
  }
  for (i in (C + 1) : M) {
    bsize_mis[i - C] ~ normal(mu_size, sd_size);
  }

  for (i in 1 : M) {
    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
                + bernoulli_logit_lpmf(y[i,] | logit_p[i,]);
    } else // s[i] == 0
    {
      for(j in 1:S){
        logit_p[i,j] = logit(mean_p) + beta * bsize_mis[i - C] + b2[j];
        target += log_sum_exp(log(p_S[j]) + bernoulli_lpmf(1 | omega)
                            // z[i] == 1
                            + bernoulli_logit_lpmf(y[i,] | logit_p[i,]),
                            bernoulli_lpmf(0 | omega) // z[i] == 0
                            //ADD categorical likelihood
                            );
      }
    }
  }
}
generated quantities {
  matrix[M, T] logit_p;
  int reach_est[M-C];
  matrix<lower=0, upper=1>[M, T] p;
  int<lower=0, upper=1> z[M];
  int<lower=C> N;

  for (i in 1 : C) {
    logit_p[i,] = to_row_vector(rep_vector(logit(mean_p) + beta * bsize[i] + b2[reach[i]], T));
  }
  for(i in (C+1):M){
    reach_est[i-C] = categorical_rng(p_S);
    logit_p[i,] = to_row_vector(rep_vector(logit(mean_p) + beta * bsize_mis[i - C] + b2[reach_est[i]],T));
  }

  p = inv_logit(logit_p);

  for (i in 1 : M) {
    if (s[i] > 0) {
      // species was detected at least once
      z[i] = 1;
    } else {
      // species was never detected
      // prob never detected given present
      real pr = prod(rep_vector(1, T) - p[i]');
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  N = sum(z);
}
