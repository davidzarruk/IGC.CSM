#' Function to invert model, so amenities, wages, productivities, and development density
#'  are chosen to match model to data.
#'
#' @param N Integer - Number of locations.
#' @param t_ij NxN matrix - Travel times across all possible locations.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param L_j Nx1 matrix - Number of workers in each location. 
#' @param Q Nx1 matrix - Floorspace prices
#' @param K Nx1 matrix - Land area
#' @param tol Int - tolerance factor
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' @param alpha Float - Utility parameter that determines preferences for
#'     consumption.
#' @param beta Float - Output elasticity wrt labor
#' @param theta Float - Commuting elasticity and migration elasticity.
#' @param lambda Float - Agglomeration force
#' @param delta Float - Decay parameter agglomeration
#' @param eta Float - Congestion force
#' @param rho Float - Decay parameter congestion
#' @param mu Float - Floorspace prod function: output elast wrt capital, 1-mu wrt land.     
#' 
#' @param zeta_intervals TODO
#' @param zeta_mult TODO
#'
#' @return Equilibrium values.
#' @export
#'
#' @examples
inversionModel = function(N,
                          L_i,
                          L_j,
                          Q,
                          K,
                          t_ij,
                          alpha=0.7,
                          beta=0.7,
                          theta=7,
                          delta=0.3585,
                          rho=0.9094,
                          lambda=0.01,
                          epsilon=0.01,
                          mu=0.3,
                          eta=0.1548,
                          nu_init=0.005,
                          tol=10^-10,
                          maxiter=100){
  
  # Formatting of input data
  L_i = array(unlist(L_i), dim(L_i))
  L_j = array(unlist(L_j), dim(L_j))
  K = array(unlist(K), dim(K))  
  Q = array(unlist(Q), dim(Q))
  t_ij = array(unlist(t_ij), dim(t_ij))  
  
  # Initialization
  w_init=array(1, dim=c(N,1))
  
  # Transformation of travel times to trade costs
  D = commuting_matrix(t_ij=t_ij, 
                       epsilon=epsilon)
  tau = D$tau
  
  # Finding the wages that match the data
  WI = wages_inversion(N=N,
                       w_init=w_init,
                       theta=theta,
                       tau=tau,
                       L_i=L_i,
                       L_j=L_j,
                       nu_init=nu_init,
                       tol=tol,
                       maxiter=maxiter)
  
  # Equilibrium wages
  w = WI$w
  w_tr = WI$w_tr
  W_i = WI$W_i
  lambda_ij_i = WI$lambda_ij_i
  
  # Average income
  Inc = av_income_simple(lambda_ij_i=lambda_ij_i,
                         w_tr = w_tr
  )
  y_bar = Inc$y_bar
  
  
  #Density of development
  DensD = density_development(Q=Q,
                              K=K,
                              w=w,
                              L_j=L_j,
                              y_bar=y_bar,
                              L_i=L_i,
                              beta=beta,
                              alpha=alpha,
                              mu=mu
  )
  Q_mean = DensD$Q_mean
  Q_norm = DensD$Q_norm
  FS_f = DensD$FS_f
  FS_r = DensD$FS_r
  FS = DensD$FS
  varphi = DensD$varphi
  ttheta = FS_f/FS
  
  #Productivities
  Prod = productivity(N=N,
                      Q=Q,
                      w=w,
                      L_j=L_j,
                      K=K,
                      t_ij = t_ij,
                      delta=delta,
                      lambda=lambda,
                      beta=beta
  )
  A = Prod$A
  a = Prod$a
  
  # Amenities
  AM = living_amenities_simple(theta=theta,
                               N=N,
                               L_i=L_i,
                               W_i=W_i,
                               Q=Q,
                               K=K,
                               alpha=alpha,
                               t_ij=t_ij,
                               rho=rho,
                               eta=eta
  )
  
  B = AM$B
  b = AM$b  
  
  # Save and export
  Q_alpha = Q_norm^(1-alpha)
  u = array_operator(array_operator(W_i,Q_alpha,'/'),B,'*')
  U = (sumDims(u,1))^(1/theta)
  
  return(list(A=A, a=a, u=u, B=B, b=b, w=w, varphi=varphi, U=U, Q_norm=Q_norm, ttheta=ttheta))
}