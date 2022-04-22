## 10. Estimation of survival, recruitment and population size using
##     the Jolly-Seber model
# 10.5. Models with individual capture heterogeneity

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Read data
## The data generation code is in bpa-code.txt, available at
## http://www.vogelwarte.ch/de/projekte/publikationen/bpa/complete-code-and-data-files-of-the-book.html
stan_data <- read_rdump("BPA/Ch.10/js_super_indran.data.R")

y2<-stan_data$y
loc<-rep(0,stan_data$M)
for(i in 1:stan_data$M){
  for(j in 2:stan_data$n_occasions){
    if(sum(y2[i,1:(j-1)]) == 2 & y2[i,j] == 1){
      y2[i,j] = 0;
      loc[i] = 1;
    }
  }
}
stan_data$y2 = y2
stan_data$loc = loc

## Parameters monitored
params <- c("sigma2", "psi", "mean_p", "mean_phi",
            "N", "Nsuper", "b", "B")

## MCMC settings
ni <- 200#30000
nt <- 1#28
nb <- 100#2000
nc <- 4

## Initial values
inits <- lapply(1:nc, function(i)
    list(mean_phi = runif(1, 0, 1),
         mean_p = runif(1, 0, 1),
         sigma = runif(1, 0, 1),
         beta = runif(stan_data$n_occasions, 0, 1)))

## Call Stan from R
js_ran <- stan("BPA/Ch.10/js_super_indran.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               seed = 2, control = list(adapt_delta = 0.9),
               open_progress = FALSE)
## lp__ of this model may have a small effective sample size.
print(js_ran, digits = 3)


## Parameters monitored
params <- c("sigma2", "psi", "p", "phi",
            "N", "Nsuper", "b", "B")
## Initial values
inits <- lapply(1:nc, function(i)
  list(phi_1 = runif(1, 0, 1),
       p_1 = runif(1, 0, 1),
       sigma = runif(1, 0, 1),
       beta = c(0,rnorm(stan_data$n_occasions-1, 0, 1))))

## Call Stan from R
js_ran <- stan("BPA/Ch.10/js_super_indran_thomas.stan",
               data = stan_data, init = inits, pars = params,
               chains = nc, iter = ni, warmup = nb, thin = nt,
               seed = 2, control = list(adapt_delta = 0.9),
               open_progress = FALSE)
## lp__ of this model may have a small effective sample size.
print(js_ran, digits = 3)
