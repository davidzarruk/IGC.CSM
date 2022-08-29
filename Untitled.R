
library(tidyverse)
library(gridExtra)

simulations = function(estimation, N, effect, beta0, beta1, beta2, interac, threshold){
  x = rnorm(n = N)
  pBO = pnorm(x)
  BO = sapply(pBO, FUN = function(x) return(rbinom(n = 1, prob = x, size = 1)))
  D = rbinom(n = N, prob = 0.50, size = 1)
  y = beta0 + effect*D + (effect+interac)*D*BO + beta1*BO + beta2*x + rnorm(n=N)
  
  mod=lm(BO ~ x, data=df)
  BOprob = mod$coefficients[1] + mod$coefficients[2]*x
  BOpred = ifelse(BOprob>threshold, 1, 0)
  
  df = data.frame(x=x, BO=BO, BOprob=BOprob, BOpred=BOpred, D=D, y=y)
  ggplot() + 
    geom_point(data = df %>% filter(BO==0), aes(x=x, y=BO), color='blue') +
    geom_point(data = df %>% filter(BO==1), aes(x=x, y=BO), color='red') +
    geom_line(data = df, aes(x=x, y=BOpred))
  
  if(estimation == 1){
    df = df %>% mutate(BOapprox = ifelse(D==1, BOpred, BO))
  } else if(estimation == 2){
    df = df %>% mutate(BOapprox = BOpred)
  } else{
    df = df %>% mutate(BOapprox = BO)
  }
  
  model = lm(y ~ D + D*BOapprox + BOapprox + x, data = df)
  
  return(model$coefficients)
}



simulations_simple = function(estimation, N, effect, beta0, beta1, beta2, interac, prob_error){
  x = rnorm(n = N)
  pBO = pnorm(x)
  BO = sapply(pBO, FUN = function(x) return(rbinom(n = 1, prob = x, size = 1)))
  D = rbinom(n = N, prob = 0.50, size = 1)
  y = beta0 + effect*D + (effect+interac)*D*BO + beta1*BO + rnorm(n=N)
  
  BOpred = rbinom(n = N, prob = prob_error, size = 1)*(1 - BO) + (1-rbinom(n = N, prob = prob_error, size = 1))*BO
  
  df = data.frame(x=x, BO=BO, BOpred=BOpred, D=D, y=y)
  
  if(estimation == 1){
    df = df %>% mutate(BOapprox = ifelse(D==1, BOpred, BO))
  } else if(estimation == 2){
    df = df %>% mutate(BOapprox = BOpred)
  } else{
    df = df %>% mutate(BOapprox = BO)
  }
  
  model = lm(y ~ D + D*BOapprox + BOapprox, data = df)
  return(model$coefficients)
}









N=10000
effect = 1
beta0 = 0.1
beta1 = 0.5
beta2 = 0.3
interac = 0.5
Nsims = 200

params_approx_05 = array(0, dim=c(Nsims, 4))
params_approx_good = array(0, dim=c(Nsims, 4))
params_approx_bad = array(0, dim=c(Nsims, 4))
params_approx_08 = array(0, dim=c(Nsims, 4))
params_true = array(0, dim=c(Nsims, 4))
for(i in 1:Nsims){
  print(i)
  # params_approx_05[i,] = simulations(estimation = TRUE, N, effect, beta0, beta1, beta2, interac, threshold=0.5)
  # params_approx_08[i,] = simulations(estimation = TRUE, N, effect, beta0, beta1, beta2, interac, threshold=0.2)
  # params_approx_good[i,] = simulations(estimation = 2, N, effect, beta0, beta1, beta2, interac, threshold=0.5)
  # params_approx_bad[i,] = simulations(estimation = 1, N, effect, beta0, beta1, beta2, interac, threshold=0.5)
  # params_true[i,] = simulations(estimation = 3, N, effect, beta0, beta1, beta2, interac, threshold=0.5)
  params_approx_good[i,] = simulations_simple(estimation = 2, N, effect, beta0, beta1, beta2, interac, prob_error=0.5)
  params_approx_bad[i,] = simulations_simple(estimation = 1, N, effect, beta0, beta1, beta2, interac, prob_error=0.5)
  params_true[i,] = simulations_simple(estimation = 3, N, effect, beta0, beta1, beta2, interac, prob_error=0.5)
}
 
g_effect_true = ggplot() + geom_histogram(aes(x = params_true[, 2])) +xlim(0.8,1.7) + geom_vline(xintercept = effect)
g_interac_true = ggplot() +  geom_histogram(aes(x = params_true[, 2]+params_true[, 3])) + xlim(0.9,3) + geom_vline(xintercept = effect+interac)

g_effect_approx_good = ggplot() +  geom_histogram(aes(x = params_approx_good[, 2])) + xlim(0.8,1.7) + geom_vline(xintercept = effect)
g_interac_approx_good = ggplot() +  geom_histogram(aes(x = params_approx_good[, 2]+params_approx_good[, 3])) + xlim(0.9,3) + geom_vline(xintercept = effect+interac)

g_effect_approx_bad = ggplot() +  geom_histogram(aes(x = params_approx_bad[, 2])) + xlim(0.8,1.7) + geom_vline(xintercept = effect)
g_interac_approx_bad = ggplot() +  geom_histogram(aes(x = params_approx_bad[, 2]+params_approx_bad[, 3])) + xlim(0.9,3) + geom_vline(xintercept = effect+interac)

g_effect_approx_05 = ggplot() +  geom_histogram(aes(x = params_approx_05[, 2])) + xlim(0.8,1.7) + geom_vline(xintercept = effect)
g_interac_approx_05 = ggplot() +  geom_histogram(aes(x = params_approx_05[, 3])) + xlim(0.2,0.6) + geom_vline(xintercept = interac)

g_effect_approx_08 = ggplot() +  geom_histogram(aes(x = params_approx_08[, 2])) + xlim(0.8,1.7) + geom_vline(xintercept = effect)
g_interac_approx_08 = ggplot() +  geom_histogram(aes(x = params_approx_08[, 3])) + xlim(0.2,0.6) + geom_vline(xintercept = interac)

grid.arrange(arrangeGrob(g_effect_true, g_interac_true, ncol=2),
             arrangeGrob(g_effect_approx_good, g_interac_approx_good, ncol=2),
             arrangeGrob(g_effect_approx_bad, g_interac_approx_bad, ncol=2))

