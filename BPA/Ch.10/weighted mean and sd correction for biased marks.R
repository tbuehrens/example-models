library(tidyverse)
library(Hmisc)
logit<-function(x){log(x/(1-x))}
ilogit<-function(x){exp(x)/(1+exp(x))}
w_mean<-function(x,w){sum((w/sum(w))*x)}
w_sd<-function(x,w,w_mean){sqrt(sum(w*(x-w_mean)^2)/(sum(w)*((length(w)-1)/length(w))))}


reps<-1000
pars<-as_tibble(matrix(NA,ncol=7))%>%
  dplyr::rename(s_mean = V1, s_sd = V2, w_mean = V3, w_sd = V4, w_sd2 = V5, mean = V6, sd = V7)
for(i in 1: reps){
  n<-10000
  x<-rnorm(n,-3,1)
  a<--1
  b<-0.8
  lo<-a+b*x
  p<-ilogit(lo)

  dat<-tibble(x,lo,p)%>%
    mutate(
      mark=rbernoulli(n,p),
      #w=(1-p)/p,
      w = (1-ilogit(b*x))/ilogit(b*x)
      )

  #true mean
  dat%>%
    summarise(mean=mean(x), sd = sd(x))
  #sample mean
  tpars<-NULL
  tpars<-dat%>%
    filter(mark==1)%>%
    summarise(s_mean=mean(x),
              s_sd = sd(x),
              w_mean = w_mean(x,w),
              w_sd = w_sd(x,w,w_mean),
              w_sd2 = sqrt(wtd.var(x,w))
              )%>%
    bind_cols(dat%>%
                summarise(mean=mean(x), sd = sd(x))
              )
  pars<-bind_rows(pars,tpars)
}

pars<-pars%>%
  filter(!is.na(sd))%>%
  mutate(err_mean = w_mean-mean,
         err_sd = w_sd - sd,
         err_s_mean = s_mean-mean,
         err_s_sd = s_sd-sd
         )

#check my weighted variance formula

plot(w_sd~w_sd2,data=pars)
abline(a=0,b=1)


hist(pars$err_mean)
hist(pars$err_sd)

quantile(pars$err_mean,c(0.025,0.25,0.5,0.75,0.975))
quantile(pars$err_sd,c(0.025,0.25,0.5,0.75,0.975))

plot(pars$w_mean~pars$mean)
plot(pars$w_sd~pars$sd)

hist(pars$err_s_mean)
hist(pars$err_s_sd)

quantile(pars$err_s_mean,c(0.025,0.25,0.5,0.75,0.975))
quantile(pars$err_s_sd,c(0.025,0.25,0.5,0.75,0.975))

plot(pars$s_mean~pars$mean)
plot(pars$s_sd~pars$sd)
