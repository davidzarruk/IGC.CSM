
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

# Parameters
alpha1 = 0.7
theta1 = 7
eta1 = 1.5
kappa1 = 2
xi1 = 1.8
nu_init = 0.005

# Number of firms
beta0 = 0.7
F = 1
sigma0 = 6

# Basic Settings
tol = 1e-6
maxiter=10
zeta = 0.0001
endo_Lr = 1
z_init=10^-4
epsilon = 0.01
mu0 = 0.3
delta0 = 0.3585
rho0 = 0.9094
eta0 = 0.1548

# Test trade costs
theta1=7

# Travel time
speed = 10

# Iceberg commuting cost
lambda = 0.01

# Workers and population
L_bar  = 1
H_bar_rest = 18

#----------------------------#
#      (2) Solve Models      #
#----------------------------#

zeta = 0.1
# Invert model
inversion_m_bl  = inversionModel(N=N,
                                L_i=L_i,
                                L_j=L_j,
                                Q=Q,
                                K=K,
                                t_ij=t_ij,
                                zeta=zeta,
                                z_init=z_init,
                                alpha=alpha1,
                                beta=beta0,
                                theta=theta1,
                                delta=delta0,
                                rho=rho0,
                                lambda=lambda,
                                epsilon=epsilon,
                                mu=mu0,
                                eta=eta0,
                                maxiter=200)


zeta = 0.0001
# Solve model
results_m_bl  = solveModel(N=N,
                          L_i=L_i,
                          L_j=L_j,
                          varphi=inversion_m_bl$varphi,
                          t_ij=t_ij,
                          K=K,
                          a=inversion_m_bl$a,
                          b=inversion_m_bl$b,
                          maxiter=10,
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
                          theta_eq=inversion_m_bl$ttheta)
