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


## Parameters monitored
params <- c("sigma2", "psi", "mean_p", "mean_phi",
            "N", "Nsuper", "b", "B")

## MCMC settings
ni <- 300#30000
nt <- 1#28
nb <- 100#2000
nc <- 4

# ## Initial values
# inits <- lapply(1:nc, function(i)
#     list(mean_phi = runif(1, 0, 1),
#          mean_p = runif(1, 0, 1),
#          sigma = runif(1, 0, 1),
#          beta = runif(stan_data$n_occasions, 0, 1)))
#
# ## Call Stan from R
# js_ran <- stan("BPA/Ch.10/js_super_indran.stan",
#                data = stan_data, init = inits, pars = params,
#                chains = nc, iter = ni, warmup = nb, thin = nt,
#                seed = 2, control = list(adapt_delta = 0.8),
#                open_progress = FALSE)
# ## lp__ of this model may have a small effective sample size.
# #print(js_ran, digits = 3)
# res<-extract(js_ran)
# quantile(res$Nsuper,c(0.025,0.25,0.5,0.75,0.975))
# apply(res$phi[,1,],2,function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975)))


y2<-stan_data$y
loc<-rep(0,stan_data$M)
loc_2=T

non_augment<-length(rowSums(stan_data$y)[rowSums(stan_data$y)>0])
loc_indexes=sort(sample(1:non_augment, round(0.5*non_augment)))

if(loc_2==T){
  for(i in 1:stan_data$M){
    if(sum(y2[i,]) == 1 & i %in% loc_indexes){
      loc[i] = 1;
    }
    for(j in 2:stan_data$n_occasions){
      if(sum(y2[i,1:j]) == 2 & y2[i,j] == 1 & i %in% loc_indexes){
        y2[i,j] = 0;
        loc[i] = 1;
      }
    }
  }
  stan_data$y = y2
  stan_data$loc = loc
}else{
  stan_data$y = stan_data$y
  stan_data$loc = rep(0,stan_data$M)
}

stan_data$ss = non_augment
stan_data$length = rnorm(non_augment,0,1)


# ## Parameters monitored
# params <- c("psi", "p", "phi",
#             "N", "Nsuper", "b", "B","chi")
# ## Initial values
# inits <- lapply(1:nc, function(i)
#   list(phi_1 = runif(1, 0, 1),
#        p_1 = runif(1, 0, 1),
#        sigma = runif(1, 0, 1),
#        beta = c(0,rnorm(stan_data$n_occasions-1, 0, 1))))
#
# ## Call Stan from R
# js_ran2 <- stan("BPA/Ch.10/js_super_indran_thomas.stan",
#                 data = stan_data, init = inits, pars = params,
#                 chains = nc, iter = ni, warmup = nb, thin = nt,
#                 seed = 2, control = list(adapt_delta = 0.8),
#                 open_progress = FALSE)
# ## lp__ of this model may have a small effective sample size.
# #print(js_ran2, digits = 3)
#
# res2<-extract(js_ran2)
# quantile(res2$Nsuper,c(0.025,0.25,0.5,0.75,0.975))
# quantile(res2$sigma2,c(0.025,0.25,0.5,0.75,0.975))
# apply(res2$phi[,1,],2,function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975)))
# apply(res2$p[,1,],2,function(x) quantile(x,c(0.025,0.25,0.5,0.75,0.975)))
# apply(res2$chi,2:3,median)



# spatial
library(gtools)
stan_data$S<-4
p_vec<-rdirichlet(1,alpha=rep(1,stan_data$S))
stan_data$x_mat<-t(rmultinom(stan_data$ss,1,p_vec))

#add in variable period length
stan_data$time = rep(1,stan_data$n_occasions-1)

## Parameters monitored
params <- c("psi", "p", "phi",
            "N", "Nsuper", "b", "B","chi")
## Initial values
inits <- lapply(1:nc, function(i)
  list(phi_1 = runif(1, 0, 1),
       p_1 = runif(1, 0, 1),
       sigma = runif(1, 0, 1),
       beta = c(0,rnorm(stan_data$n_occasions-1, 0, 1))))

## Call Stan from R
js_ran2 <- stan("BPA/Ch.10/js_super_indran_thomas_test_v2_space.stan",
                data = stan_data, init = inits, pars = params,
                chains = nc, iter = ni, warmup = nb, thin = nt,
                seed = 2, control = list(adapt_delta = 0.8),
                open_progress = FALSE)
## lp__ of this model may have a small effective sample size.
#print(js_ran2, digits = 3)
res2<-extract(js_ran2)
quantile(res2$Nsuper,c(0.025,0.25,0.5,0.75,0.975))
write.csv(summary(js_ran2)$summary,"BPA/Ch.10/summary.csv")

