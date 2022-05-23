
#-----------------#
#   Packages      #
#-----------------#

install.packages('devtools')
library('devtools')

install_github("davidzarruk/IGCities", force = TRUE)
library(IGCities)

#-----------------#
#   Parameters    #
#-----------------#

N = 875
S = 3

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

# Test trade costs
theta1=7

# Travel time
speed = 10

# Iceberg commuting cost
lambda = 0.01

# Workers and population
L_bar  = 1
H_bar_rest = 18

# Loaded parameters
H_bar = read.csv("data/H_bar.csv")
tau = read.csv("data/tau.csv")
L_j_data = read.csv("data/L_j_data.csv")
lambda_is_i = read.csv("data/lambda_is_i.csv")
lambda_i = read.csv("data/lambda_i.csv")


#----------------------------#
#      (2) Solve Models      #
#----------------------------#

zeta = 0.1
# Invert model
inversion_m_bl  = inversionModel_Eff(N=N,
                                     S=S,
                                     L_bar=L_bar,
                                     H_bar=H_bar,
                                     H_bar_rest=H_bar_rest,
                                     tau=tau,
                                     lambda_i=lambda_i,
                                     lambda_is_i=lambda_is_i,
                                     L_j_data=L_j_data,
                                     zeta=zeta,
                                     z_init=z_init,
                                     tol=tol,
                                     maxiter=maxiter,
                                     alpha1=alpha1,
                                     beta0=beta0,
                                     theta1=theta1,
                                     eta1=eta1,
                                     kappa1=kappa1,
                                     sigma0=sigma0,
                                     xi1=xi1,
                                     F=F);


zeta = 0.0001
# Solve model
results_m_bl  = solveModel1_Eff(N=N,
                                S=S,
                                L_bar=L_bar,
                                H_bar=H_bar,
                                H_bar_rest=H_bar_rest,
                                tau=tau,
                                A=inversion_m_bl$A,
                                u_eq=inversion_m_bl$u,
                                B=inversion_m_bl$B,
                                w_eq=inversion_m_bl$w,
                                lambda_i=lambda_i,
                                lambda_is_i=lambda_is_i,
                                zeta=zeta,
                                tol=tol,
                                maxiter=maxiter,
                                endo_Lr=endo_Lr,
                                alpha1=alpha1,
                                beta0=beta0,
                                theta1=theta1,
                                eta1=eta1,
                                kappa1=kappa1,
                                xi1=xi1,
                                sigma0=sigma0,
                                F=F)
