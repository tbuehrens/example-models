// JS model using the superpopulation parameterization;


//1 run random walk and OG version to verify same results (DONE)
//2 modify data to include LOC and repeat step 1 (DONE)
//3 use Dan email to check calcs for B* instead of B (including fish born + died before available for capture within a period)
//4 flow???
functions {
  // These functions are derived from Section 12.3 of
  // Stan Modeling Language User's Guide and Reference Manual

  /**
   * Return a integer value of first capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of first capture occasion
   */
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i)){
      if (y_i[k]){
        return k;
      }
    }
    return 0;
  }

  /**
   * Return a integer value of last capture occasion
   *
   * @param y_i Integer array of capture history
   *
   * @return Integer value of last capture occasion
   */
  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      int k;
      k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }

  /**
   * Return a matrix of uncaptured probabilities
   *
   * @param p   Matrix of detection probabilities for each individual
   *            and capture occasion
   * @param phi Matrix of survival probabilities for each individual
   *            and capture occasion
   *
   * @return Matrix of uncaptured probabilities
   */
  matrix prob_uncaptured(matrix p, matrix phi, int[] loc, int[] last) {
    int n_ind = rows(p);
    int n_occasions = cols(p);
    matrix[n_ind, n_occasions] chi;

    for (i in 1 : n_ind) {
      chi[i, n_occasions] = 1.0;
      for (t in 1 : (n_occasions - 1)) {
        int t_curr = n_occasions - t;
        int t_next = t_curr + 1;
        chi[i, t_curr] = (1 - phi[i, t_curr]) + phi[i, t_curr] * (1 - p[i, t_next]) * chi[i, t_next];
      }
      for (t in 1: (n_occasions - 1)) {
        int t_curr = n_occasions - t;
        if(loc[i] == 1){
          if(t_curr >= last[i]){
            chi[i, t_curr] = 1; //set prob uncaptured to 1.0 after loss on capture
          }
        }
      }
    }
    return chi;
  }

  row_vector prob_uncaptured_ind(row_vector p, row_vector phi, int loc, int last) {
    int n_occasions = cols(p);
    row_vector[n_occasions] chi;
    chi[n_occasions] = 1.0;
    for (t in 1 : (n_occasions - 1)) {
      int t_curr = n_occasions - t;
      int t_next = t_curr + 1;
      chi[t_curr] = (1 - phi[t_curr]) + phi[t_curr] * (1 - p[t_next]) * chi[t_next];
    }
    for (t in 1: (n_occasions - 1)) {
      int t_curr = n_occasions - t;
      if(loc == 1){
        if(t_curr >= last){
          chi[t_curr] = 1; //set prob uncaptured to 1.0 after loss on capture
        }
      }
    }
    return chi;
  }

  /**
   * Calculate log likelihood of a Jolly-Seber model
   * under the superpopulation parameterization
   *
   * @param y     Integer array of capture history
   * @param first Integer array of first capture occasions
   * @param last  Integer array of last capture occasions
   * @param p     Matrix of detection probabilities
   * @param phi   Matrix of survival probabilities
   * @param psi   Real value of inclusion probability
   * @param nu    Vector of entry probabilities
   * @param chi   Matrix of uncapture probabilities
   */
  void js_super_lp(int[,] y,
                   int[] first,
                   int[] last,
                   //matrix p,
                   //matrix phi,
                   int[] location,
                   vector length,
                   vector length_est,
                   int S,
                   vector time,
                   int[] loc,

                   vector spatial_vec,

                   real p_1,
                   real b1_p,
                   vector b2_p,
                   vector eps,
                   real sigma,
                   vector eps_p_t,
                   real sigma_p_t,

                   real phi_1,
                   real b1_phi,
                   vector b2_phi,
                   vector  eps_phi_t,
                   real sigma_phi_t,

                   real psi,
                   vector nu
                   ) {
    int n_ind = dims(y)[1];
    int n_occasions = dims(y)[2];
    vector[n_occasions] qnu = 1.0 - nu;

    int ss = dims(location)[1];
    matrix[n_ind, n_occasions] p; //prob capture
    matrix[n_ind, n_occasions-1] phi; //prob survival
    for (i in 1 : n_ind) {
      vector[n_occasions] qp;
      row_vector[n_occasions] chi;//prob non-recovered

      if (first[i]) {
        // Observed
        p[i,1] = inv_logit(logit(p_1) + b1_p * length[i] + b2_p[location[i]] + eps[i] * sigma);
        for (t in 2 : (n_occasions)){
          p[i,t] = inv_logit(logit(p[i,t-1]) + eps_p_t[t-1] * sigma_p_t);
        }

        phi[i,1] = inv_logit(logit(phi_1) + b1_phi * length[i] + b2_phi[location[i]]) ^ time[1];
        for (t in 2 : (n_occasions-1)){
          phi[i,t] = inv_logit(logit(phi[i,t-1]^time[t-1]^-1) + eps_phi_t[t-1] * sigma_phi_t) ^ time[t];
        }

        qp = 1.0 - p[i]';
        // Uncaptured probability
        chi = prob_uncaptured_ind(p[i], phi[i], loc[i], last[i]);

        // Included
        1 ~ bernoulli(psi);

        // Until first capture
        if (first[i] == 1) {
          1 ~ bernoulli(nu[1] * p[i, 1]); //may need to add in prob of survival BEFORE first capture or maybe use equation for births
        } else {
          // first[i] >= 2
          vector[first[i]] lp;

          // Entered at 1st occasion
          lp[1] = bernoulli_lpmf(1 | nu[1])
                  + bernoulli_lpmf(1 | prod(qp[1 : (first[i] - 1)]))
                  + bernoulli_lpmf(1 | prod(phi[i, 1 : (first[i] - 1)]))
                  + bernoulli_lpmf(1 | p[i, first[i]]);
          // Entered at t-th occasion (1 < t < first[i])
          for (t in 2 : (first[i] - 1)) {
            lp[t] = bernoulli_lpmf(1 | prod(qnu[1 : (t - 1)]))
                    + bernoulli_lpmf(1 | nu[t])
                    + bernoulli_lpmf(1 | prod(qp[t : (first[i] - 1)]))
                    + bernoulli_lpmf(1 | prod(phi[i, t : (first[i] - 1)]))
                    + bernoulli_lpmf(1 | p[i, first[i]]);
          }
          lp[first[i]] = bernoulli_lpmf(1 | prod(qnu[1 : (first[i] - 1)]))
                         + bernoulli_lpmf(1 | nu[first[i]])
                         + bernoulli_lpmf(1 | p[i, first[i]]);
          target += log_sum_exp(lp);
        }
        // Until last capture
        for (t in (first[i] + 1) : last[i]) {
          1 ~ bernoulli(phi[i, t - 1]); // Survived
          y[i, t] ~ bernoulli(p[i, t]); // Capture/Non-capture
        }
        // Subsequent occasions
        1 ~ bernoulli(chi[last[i]]);
      } else {
        vector[S] lps = log(spatial_vec);
        for(s in 1:S){
        // Never observed
          vector[n_occasions + 1] lp;
          p[i,1] = inv_logit(logit(p_1) + b1_p * length_est[i-ss] + b2_p[s] + eps[i] * sigma);
          for (t in 2 : (n_occasions)){
            p[i,t] = inv_logit(logit(p[i,t-1]) + eps_p_t[t-1] * sigma_p_t);
          }
          phi[i,1] = inv_logit(logit(phi_1) + b1_phi * length_est[i-ss] + b2_phi[s]) ^ time[1];
          for (t in 2 : (n_occasions-1)){
            phi[i,t] = inv_logit(logit(phi[i,t-1]^time[t-1]^-1) + eps_phi_t[t-1] * sigma_phi_t) ^ time[t];
          }

          qp = 1.0 - p[i]';
          // Uncaptured probability
          chi = prob_uncaptured_ind(p[i], phi[i], loc[i], last[i]);

          // Entered at 1st occasion, but never captured
          lp[1] = bernoulli_lpmf(1 | psi) + bernoulli_lpmf(1 | nu[1])
                  + bernoulli_lpmf(0 | p[i, 1]) + bernoulli_lpmf(1 | chi[1]);
          // Entered at t-th occation (t > 1), but never captured
          for (t in 2 : n_occasions) {
            lp[t] = bernoulli_lpmf(1 | psi)
                    + bernoulli_lpmf(1 | prod(qnu[1 : (t - 1)]))
                    + bernoulli_lpmf(1 | nu[t]) + bernoulli_lpmf(0 | p[i, t])
                    + bernoulli_lpmf(1 | chi[t]);
          }
          // Never captured
          lp[n_occasions + 1] = bernoulli_lpmf(0 | psi);
          lps[s] += log_sum_exp(lp);
        }
        target += log_sum_exp(lps);
      }
    }
  }
}
data {
  int<lower=0> ss; //un-augmented sample size (top rows of capture history matrix above all zero augmented rows)
  vector[ss] length; //z_scored lengths for observed capture histories
  int S; //count of spatial sections
  int location [ss];//design matrix for which section an observed fish was caught in
  int<lower=0> M; // Augmented sample size
  int<lower=0> n_occasions; // Number of capture occasions
  vector<lower=0>[n_occasions-1] time;//period length (assumed same for all sectons currently); no length for first period (before first event)
  int loc[M];//if loc == 1, fish lost on capture
  int<lower=0, upper=1> y[M, n_occasions]; // Augmented capture-history
}
transformed data {
  int x_vec[S];//counts of observed fish by locatiom
  int<lower=0, upper=n_occasions> first[M];
  int<lower=0, upper=n_occasions> last[M];
  for(i in 1:S){
    x_vec[i] = 0;
  }
  for(i in 1:ss){
    x_vec[location[i]] += 1;
  }
  for (i in 1 : M) {
    first[i] = first_capture(y[i]);
  }
  for (i in 1 : M) {
    last[i] = last_capture(y[i]);
  }
}
parameters {
  //imputed covariates
  vector[M-ss] z_eps_length; //z-scored weight residuals for unobserved fish
  real mean_length; //mean length of all fish
  real<lower=0> sd_length; //sd length of all fish
  simplex[S] spatial_vec;//proportions of fish by location for all fish
  //survival
  real<lower=0, upper=1> phi_1; // survival first period
  vector[n_occasions-2] eps_phi_t; //temporal process errors in phi
  real<lower=0> sigma_phi_t; //temporal process error SD for phi
  vector[S-1] eps_phi_s; //spatial process errors in phi
  real<lower=0> sigma_phi_s; //spatial process error SD for phi
  real b1_phi; //coefficient for effect of length on prob of survival
  //capture
  real<lower=0, upper=1> p_1; // capture first period
  vector[M] eps;//individual random effect prob of capture
  real<lower=0> sigma;//individual random effect SD prob of capture
  vector[n_occasions-1] eps_p_t; //temporal process errors in p
  real<lower=0> sigma_p_t; //temporal process error SD for p
  real b1_p; //coefficient for effect of length on prob of capture
  vector[S-1] eps_p_s; //spatial process errors in p
  real<lower=0> sigma_p_s; //spatial process error SD for p
  //inclusion
  real<lower=0, upper=1> psi; // Inclusion probability
  //entry
  vector[n_occasions-1] eps_b_t; //temporal process errors in births
  real<lower=0> sigma_b_t; //temporal process error SD for births
}
transformed parameters {
  //imputed covariates
  vector[M-ss] length_est; //lengths of uncaptured
  //capture
  vector[S] b2_p;//section-specific survival adjustment
  //surviva/;
  vector[S] b2_phi;//section-specific survival adjustment
  //entry
  vector[n_occasions] beta; //additive log ratio component for prob of entry
  simplex[n_occasions] b; // Entry probability
  vector<lower=0, upper=1>[n_occasions] nu;//conditional prob entry
  //other

  length_est = mean_length + z_eps_length * sd_length;

  b2_phi[1] = 0;
  b2_p[1] = 0;
  for(s in 2:S){
    b2_phi[s] = b2_phi[s-1] + eps_phi_s[s-1] * sigma_phi_s;
    b2_p[s] = b2_p[s-1] + eps_p_s[s-1] * sigma_p_s;
  }

  // Additive log ratio random walk prior for entry probabilities w
  beta[1] = 0; //no need for time offset for first period
  for (t in 2 : n_occasions) {
    beta[t] = beta[1] + sum(eps_b_t[1:(t-1)] * sigma_b_t + log(time[t-1]));
  }
  b = softmax(beta[1:n_occasions]);

  // Convert entry probs to conditional entry probs (conditional that you didn't enter prior to current period)
  {
    real cum_b = b[1];

    nu[1] = b[1];
    for (t in 2 : (n_occasions - 1)) {
      nu[t] = b[t] / (1.0 - cum_b);
      cum_b = cum_b + b[t];
    }
    nu[n_occasions] = 1.0;
  }
}
model {
  // Priors
  //covariate imputation
  z_eps_length ~ std_normal();
  mean_length ~ std_normal();
  sd_length ~ std_normal();
  //capture
  //p_1 ~ uniform(0,1);
  sigma ~ std_normal();
  eps ~ std_normal();
  sigma_p_t ~ std_normal();
  eps_p_t ~ std_normal();
  b1_p ~std_normal();
  eps_p_s ~ std_normal();
  sigma_p_s ~ std_normal();
  //survival
  //phi_1 ~ uniform(0,1);
  b1_phi ~std_normal();
  sigma_phi_t ~ std_normal();
  eps_phi_t ~ std_normal();
  sigma_phi_s ~ std_normal();
  eps_phi_s ~ std_normal();
  //entry
  sigma_b_t ~ std_normal();
  eps_b_t ~ std_normal();

  // Likelihoods
  js_super_lp(y,
              first,
              last,
              location,
              length,
              length_est,
              S,
              time,
              loc,
              spatial_vec,
              p_1,
              b1_p,
              b2_p,
              eps,
              sigma,
              eps_p_t,
              sigma_p_t,
              phi_1,
              b1_phi,
              b2_phi,
              eps_phi_t,
              sigma_phi_t,
              psi,
              nu)
  length ~ normal(mean_length,sd_length);
  x_vec ~ multinomial(spatial_vec);
}
generated quantities {
  int location_aug[M];//augmented design matrix for which section an observed fish was caught in
  matrix[M, n_occasions] p; //prob capture
  matrix[M, n_occasions-1] phi; //prob capture
  vector[M] length_aug; //lengths of whole pop
  real<lower=0> sigma2;
  int<lower=0> Nsuper; // Superpopulation size
  int<lower=0> N[n_occasions]; // Actual population size
  int<lower=0> B[n_occasions]; // Number of entries
  int<lower=0, upper=1> z[M,n_occasions]; // Deflated latent state

  sigma2 = square(sigma);

  //merge augemented and non-augmented lengths
  length_aug[1:ss] = length;
  length_aug[(ss+1):M] = length_est;


  //generate augmented x_mat, p and phi
  for(i in 1:M){
    if(i <= ss){

    }else{
       //to_row_vector(spatial_vec);
    }
  }
  for(i in 1:M){
    if(i <= ss){
      location_aug[i] = location[i];
    }else{
      location_aug[i] =  categorical_rng(spatial_vec);
    }
    p[i,1] = inv_logit(logit(p_1) + b1_p * length_aug[i] + b2_p[location_aug[i]]);
    for (t in 2 : (n_occasions)){
      p[i,t] = inv_logit(logit(p[i,t-1]) + eps_p_t[t-1] * sigma_p_t);
    }
    phi[i,1] = inv_logit(logit(phi_1) + b1_phi * length_aug[i] + b2_phi[location_aug[i]]) ^ time[1];
    for (t in 2 : (n_occasions-1)){
      phi[i,t] = inv_logit(logit(phi[i,t-1]^time[t-1]^-1) + eps_phi_t[t-1] * sigma_phi_t) ^ time[t];
    }
  }
  // Generate w[] and z[]
  for (i in 1 : M) {
    int q = 1;
    if (bernoulli_rng(psi)) {
      // Included
      z[i, 1] = bernoulli_rng(nu[1]);
      for (t in 2 : n_occasions) {
        q = q * (1 - z[i, t - 1]);
        z[i, t] = bernoulli_rng(z[i, t - 1] * phi[i, t - 1] + q * nu[t]);
      }
    } else {
      // Not included
      z[i,  : ] = rep_array(0, n_occasions);
    }
  }

  // Calculate derived population parameters
  {
    int recruit[M, n_occasions] = rep_array(0, M, n_occasions);
    int Nind[M];
    int Nalive[M];

    for (i in 1 : M) {
      int f = first_capture(z[i,  : ]);

      if (f > 0) {
        recruit[i, f] = 1;
      }
    }
    for (t in 1 : n_occasions) {
      N[t] = sum(z[ : , t]);
      B[t] = sum(recruit[ : , t]);
      //Add B* using dan email calc here...may need to loop through periods AND individuals to get this if individual heterogeneous survival probs within period
    }
    for (i in 1 : M) {
      Nind[i] = sum(z[i]);
      Nalive[i] = Nind[i] > 0;
    }
    Nsuper = sum(Nalive);// see if Nsuper2 = sum (B[1:n_occasions]) is same!
    //calc new Nsuper as sum of Bstars
  }
}
