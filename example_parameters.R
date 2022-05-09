
library(roxygen2)
setwd('/home/david/Documents/IGCities/')
rm(list = ls())

#-----------------#
#   Load data     #
#-----------------#

Data = read.csv('data/data_shriaka.csv', sep='\t')
Jobs = read.csv('data/employment_shriaka.csv', sep='\t')
Resi = read.csv('data/residents_shriaka.csv', sep='\t')

#-----------------#
#   Parameters    #
#-----------------#

N = 875
S = 3

# Floor space supply
delta0 = 0.65

# Parameters
alpha1 = 0.7
theta1 = 4
eta1 = 1.5
kappa1 = 2
xi1 = 1.8

# Number of firms
beta0 = 0.7
F = 1
sigma0 = 6

# Basic Settings
tol = 1e-10
zeta = 0.0001
endo_Lr = 1

# Test trade costs
theta1=7

# Travel time
speed = 10

# Iceberg commuting cost
lambda = 0.01

# Workers and population
L_bar  = 1

area = Data$Area/1000000
FAR = Data$far
H_bar_rest = 18
H_bar = FAR*area


#-----------------#
# Commuting costs #
#-----------------#

dist_mat = read.table("data/distance_matrix.txt", sep = '\t')
dist_mat = dist_mat[,2:876]
dist_mat = array(unlist(dist_mat), dim(dist_mat))

# Calculate travel time. Assuming speed of 15 km/h
time_mat = (dist_mat/speed)*60

# Iceberg commuting cost
tau = exp(lambda*time_mat)

#------------------------#
# Workers and population #
#------------------------#

# Workers and population
pop_t = Data$pop
workers_t = Jobs$total_empl_shriaka

workers_s = Jobs[,3:5]*L_bar
workers_s = array(unlist(workers_s), dim(workers_s))

pop = pop_t/sum(pop_t)
workers = workers_t/sum(workers_t)
lambda_is_i = Resi[,2:4]
lambda_is_i = array(unlist(lambda_is_i), dim(lambda_is_i))

res_agg_s = array(0, dim=c(1,3))
wor_agg_s = array(0, dim=c(1,3))
for(irow in 1:875){
  for(is in 1:3){
    res_agg_s[1, is] = res_agg_s[1, is] + (kronecker(pop, array(1, dim=c(1,3)))*lambda_is_i)[irow, is]
    wor_agg_s[1, is] = wor_agg_s[1, is] + (workers_s)[irow, is]
  }
}

L_j_data = (workers_s*kronecker(res_agg_s, array(1, dim=c(875,1))))/kronecker(wor_agg_s, array(1, dim=c(875,1)))

# Residents in each location
lambda_i = pop


#------------------------#
# Export example files   #
#------------------------#

write.table(tau, "data/tau.txt", sep = '\t')
write.csv(L_j_data, "data/L_j_data.csv", row.names = FALSE)
write.csv(lambda_is_i, "data/lambda_is_i.csv", row.names = FALSE)
write.csv(lambda_i, "data/lambda_i.csv", row.names = FALSE)
write.csv(H_bar, "data/H_bar.csv", row.names = FALSE)
