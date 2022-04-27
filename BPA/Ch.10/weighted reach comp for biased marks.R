library(tidyverse)
library(Hmisc)
library(gtools)
logit<-function(x){log(x/(1-x))}
ilogit<-function(x){exp(x)/(1+exp(x))}
logsumexp <- function (x) {
  y = max(x)
  y + log(sum(exp(x - y)))
}
softmax <- function (x) {
  exp(x - logsumexp(x))
}

S<-4
sigma_phi_s<-0.5
p_space<-rdirichlet(1,alpha=rep(1,S))
eps_phi_s = rnorm(S,0,1)
b_space_phi<-rep(NA,S)
b_space_phi[1] = 0;
for(s in 2:S){
  b_space_phi[s] = b_space_phi[s-1] + eps_phi_s[s-1] * sigma_phi_s;
}



p_space_all = softmax(logit(p_space) - b_space_phi)


