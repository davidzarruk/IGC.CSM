
#-----------------#
#   Packages      #
#-----------------#

install.packages('devtools')
library('devtools')

install_github("davidzarruk/IGCities", force = TRUE)
library(IGCities)

#setwd("C:/Users/romandzarate/GitHub/IGCities")
setwd('/Users/zarruk/Documents/IGCities/')

#-----------------#
#   Parameters    #
#-----------------#

rm(list = ls())

# Data
data_locations = read.csv("data/Data for model/Chars.csv")
data_times = read.csv("data/Data for model/MatrixTravelTimes_mins.csv")

L_j = as.data.frame(data_locations$t_w_vodacom) # cantidad de trabajadores en cada location j
L_i = as.data.frame(data_locations$t_r_vodacom) # cantidad de habitantes en cada location j
L_i = L_i*sum(L_j)/sum(L_i)
K = as.data.frame(data_locations$SAL_area) # tamaño del lugar
Q = as.data.frame(data_locations$price_m2)
N = dim(L_i)[1]
t_ij = as.matrix(data_times[,2:(N+1)], dim=c(N,N))
t_ij[864,] = t_ij[863,]
t_ij[,864] = t_ij[,863]
t_ij[,714] = t_ij[,713]
t_ij[714,] = t_ij[713,]
t_ij[165,] = t_ij[164,]
t_ij[,165] = t_ij[,164]
t_ij[,138] = t_ij[,137]
t_ij[138,] = t_ij[137,]

#----------------------------#
#      (2) Solve Models      #
#----------------------------#

zeta = 0.3
# Invert model
inversion_m_bl  = inversionModel(N=N,
                                L_i=L_i,
                                L_j=L_j,
                                Q=Q,
                                K=K,
                                t_ij=t_ij)

increase_prod = L_j$`data_locations$t_w_vodacom` > 7.837954e+03
a_increase = increase_prod*inversion_m_bl$a*1.1 + (1-increase_prod)*inversion_m_bl$a

zeta = 0.95
# Solve model
results_m_bl  = solveModel(N=N,
                          L_i=L_i,
                          L_j=L_j,
                          varphi=inversion_m_bl$varphi,
                          t_ij=t_ij,
                          K=K,
                          a=a_increase,
                          b=inversion_m_bl$b,
                          maxiter=500,
                          alpha=alpha1,
                          beta=beta0,
                          theta=theta1,
                          mu=mu0,
                          delta=delta0,
                          lambda=lambda,
                          rho=rho0,
                          eta=eta0,
                          epsilon=epsilon,
                          w_eq=inversion_m_bl$w,
                          u_eq=inversion_m_bl$u,
                          Q_eq=inversion_m_bl$Q_norm,
                          ttheta_eq=inversion_m_bl$ttheta)

# Shock
a = inversion_m_bl$a
p90=quantile(a, 0.9)
a_new = a*(1+0.1*(a>p90))

# Solve model after shock
results_m_bl  = solveModel(N=N,
                           L_i=L_i,
                           L_j=L_j,
                           varphi=inversion_m_bl$varphi,
                           t_ij=t_ij,
                           K=K,
                           a=a_new,
                           b=inversion_m_bl$b,
                           maxiter=500,
                           alpha=alpha1,
                           beta=beta0,
                           theta=theta1,
                           mu=mu0,
                           delta=delta0,
                           lambda=lambda,
                           rho=rho0,
                           eta=eta0,
                           epsilon=epsilon,
                           w_eq=inversion_m_bl$w,
                           u_eq=inversion_m_bl$u,
                           Q_eq=inversion_m_bl$Q_norm,
                           ttheta_eq=inversion_m_bl$ttheta)

