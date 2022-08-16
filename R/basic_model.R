#' Function to transform travel times into iceberg commuting costs
#' 
#' @param t_ij NxN matrix - Travel time matrix across locations
#' @param delta Float - Parameter that transforms travel times to commuting costs
#' 
#' @return A NxN matrix of commuting costs
#' 
#' @export
#' 
#' @examples 
commuting_matrix = function(t_ij,
                            epsilon){
  tau = exp(epsilon*t_ij)
  return(list(tau=tau))
}


#' Function to compute equilibrium wages that make the model labor in every
#' location in equal to the observed data. It finds the w's
#' such that equation (3.2) holds.
#'
#' @param N Integer - Number of locations.
#' @param w_init Initial vector of wages.
#' @param theta Float - Commuting elasticity.
#' @param tau NxN matrix - Commuting cost matrix across all locations.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param L_j Nx1 matrix - Number of workers in each location.
#' @param nu_init Float - Convergence parameter to update wages.
#'     Default nu=0.01.
#' @param tol Float - Maximum tolerable error for estimating total labor.
#'     Default tol=10^-10.
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=10000.
#'
#' @return A list with equilibrium wages and probability of workers in each
#'     location working in every other location.
#' @export
#'
#' @examples
wages_inversion = function(N,
                          w_init,
                          theta,
                          tau,
                          L_i,
                          L_j,
                          nu_init=0.05,
                          tol=10^-10,
                          maxiter=10000,
                          nu_mult=c(),
                          nu_intervals=c()){

  # Settings
  outerdiff = Inf
  w = w_init
  iter = 0
  nu = nu_init

  while(outerdiff>tol & iter<maxiter){
    # 1) Labor supply
    # Indirect utility
    w_tr = aperm(array(w, dim=c(N,1)), c(2,1));
    rep_w_tr = kronecker(w_tr^theta, array(1, dim=c(N, 1)));
    # Constructing emp` loyment shares
    w_tr_tau = array_operator(w_tr^theta, tau^(-theta), '*');
    lambda_ij_i = array_operator(w_tr_tau, sumDims2(w_tr_tau,2), '/');
    W_i = (sumDims2(w_tr_tau,2))^(1/theta);
    
    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    L_ij = array_operator(L_i, lambda_ij_i, '*')
    L_j_tr = sumDims2(L_ij, 1)
#    L_j_model = aperm(L_j_tr, c(2, 1));
    Ratio_supply = array_operator(L_j_tr, w, "/");
    w_prime = array_operator(L_j, Ratio_supply, "/");

    z_L = array_operator(w, w_prime, '-');
    w = array_operator(w*(1-nu), w_prime*nu, '+');
    w_mean = exp(mean(log(w)))
    w = w/w_mean;
    outerdiff = max(abs(z_L))
    
    iter = iter+1;

    if(iter>100){
      nu = nu_init*10;
    } else if(iter>1000){
      nu = nu_init*10;
    } else if(iter>1000000){
      break
    }

    # if(length(nu_mult) > 0){
    #   if(length(nu_mult) != length(nu_intervals)){
    #     stop(
    #       "The array nu_mult must have the same size as the array nu_intervals."
    #     )
    #   } else{
    #     nu_mult2 = c(1, nu_mult)
    #     iiter = 1
    #     while(iter > nu_intervals[iiter] & iiter <= length(nu_intervals)){
    #       iiter = iiter+1
    #     }
    #     nu=nu_init*nu_mult2[iiter]
    #   }
    # }
    print(paste(outerdiff, iter))
  }

  return(list(w=w, w_tr=w_tr, W_i=W_i, lambda_ij_i=lambda_ij_i, L_j_model=L_j_model))
}


#' Computes average income in each location, which is the weighted average of
#' the income of the people living in the location.
#'
#' @param lambda_ij_i NxN matrix - Probability of individuals in each
#'     location of working in each location.
#' @param w NxS - Wages in each location in each sector.
#'
#' @return
#' @export
#'
#' @examples
av_income = function(lambda_ij_i,
                     w_tr){

  y_bar = (sumDims2(array_operator(lambda_ij_i, w_tr, '*'), 2));
  return(list(y_bar=y_bar))
}

#' Computes residential and commercial floorspace supply and equilibrium prices.
#'
#' @param Q Nx1 array - Floorspaces prices.
#' @param K Nx1 array - Land supply.
#' @param w NxS - Wages in each location in each sector.
#' @param L_j Nx1 matrix - Number of workers in each location.
#' @param y_bar - Average income in each location.
#' @param L_i Nx1 matrix - Number of residents in each location.
#' @param beta Float - Cobb-Douglas parameter output elasticity wrt labor.
#' @param alpha Float - Utility parameter that determines preferences for
#'     consumption.
#' @param mu Float - Floorspace prod function: output elast wrt capita, 1-mu wrt land.     
#'
#' @return
#' @export
#'
#' @examples
density_development = function(Q,
                               K,
                               w,
                               L_j,
                               y_bar,
                               L_i,
                               beta,
                               alpha,
                               mu){

  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  FS_f = ((1-beta)/beta)*(array_operator(array_operator(w, L_j, '*'), Q_norm, '/'));
  FS_r = (1-alpha)*(array_operator(array_operator(y_bar, L_i, '*'),Q_norm,'/'));
  FS = FS_f+FS_r;
  varphi = array_operator(FS, K^(1-mu), '/'); 
  return(list(Q_mean = Q_mean, Q_norm = Q_norm, FS_f = FS_f, FS_r = FS_r, FS = FS, varphi=varphi))
}

#' Computes productivity levels in each location
#'
#' @param N Float - Number of locations.
#' @param w Nx1 matrix - wages in each location.
#' @param Q Nx1 matrix - Floorspace prices in each location.
#' @param L_j Nx1 matrix - Employment in each location.
#' @param K Nx1 matrix - Land in each location.
#' @param t_ij NxN matrix - Travel times matrix.
#' @param delta Float - decay parameter agglomeration.
#' @param lambda Float - agglomeration force.
#' @return
#' @export
#'
#' @examples
productivity = function(N,
                        Q,
                        w,
                        L_j,
                        K,
                        t_ij,
                        delta,
                        lambda){

  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  beta_tilde = ((1-beta)^(1-beta))*(beta^beta); 
  A = (1/beta_tilde)*(array_operator(Q_norm^(1-beta), w^beta, '*'));
  L_j_dens = (array_operator(L_j, K, '/'));
  L_j_dens_per = aperm(array(L_j_dens, dim=c(N,1)), c(2,1));
  L_j_dens_rep = kronecker(L_j_dens_per, array(1, dim=c(N, 1)));
  Upsilon = sumDims2(array_operator(exp(-delta*t_ij), L_j_dens_rep, '*'), 2);
  a = array_operator(A, Upsilon^lambda, "/");
  a = a*(L_j>0)
  return(list(A = A, a = a))
}

#' Function to estimate amenity parameters of locations where users live.
#'
#' @param theta Float - Parameter that governs the reallocation of workers across
#'     locations in the city. This parameter measures how sensible are migration
#'     flows within the city to changes in real income.
#' @param alpha Float - Para     
#' @param Q Nx1 matrix - Floor space prices.
#' @param L_i Nx1 matrix - Total residents.
#' @param W_i Nx1 matrix - Market access measure in each location.
#' @param t_ij NxN matrix - Travel times across locations.
#' @param rho Float - decay parameter for amenities.
#' @param eta Float - congestion force
#'
#' @return Matrix with the amenity distribution of living in each location.
#' @export
#'
#' @examples
living_amenities = function(theta,
                            N,
                            L_i,
                            W_i,
                            Q,
                            alpha,
                            t_ij,
                            rho,
                            eta){
  Q_mean = exp(mean(log(Q)));
  Q_norm = Q/Q_mean;
  L_i_mean = exp(mean(log(L_i)));
  L_i_norm = L_i/L_mean;
  W_i_mean = exp(mean(log(W_i)));
  W_i_norm = W_i/W_i_mean;
  B = array_operator(array_operator(L_i_norm^(1/theta), Q_norm^(1-alpha), '*'), W_i_norm^((-1)/theta), '*');
  L_i_dens = (array_operator(L_i, K, '/'));
  L_i_dens_per = aperm(array(L_i_dens, dim=c(N,1)), c(2,1));
  L_i_dens_rep = kronecker(L_i_dens_per, array(1, dim=c(N, 1)));
  Omega = sumDims2(array_operator(exp(-rho*t_ij), L_i_dens_rep, '*'), 2);  
  b = array_operator(B, Omega^(-eta), "/");
  b = b*(L_i>0);
  return(list(B = B, b = b))
}


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
                          zeta,
                          z_init,
                          tol=10^-10,
                          maxiter=1000,
                          alpha,
                          beta,
                          theta,
                          delta,
                          rho,
                          lambda,
                          epsilon,
                          mu,
                          nu_init=0.005,
                          nu_intervals = c(),
                          nu_mult = c(),
                          zeta_intervals = c(),
                          zeta_mult = c()){

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
                     maxiter=maxiter,
                     nu_intervals = nu_intervals,
                     nu_mult = nu_mult)

  # Equilibrium wages
  w = WI$w
  w_tr = WI$w_tr
  W_i = WI$W_i
  lambda_ij_i = WI$lambda_ij_i
  
  # Average income
  Inc = av_income(lambda_ij_i=lambda_ij_i,
                   w_tr = w_tr
                  )
  ybar = Inc$y_bar
  
  
  #Density of development
  DensD = density_development(Q=Q,
                              K=K,
                              w=w,
                              L_j=L_j,
                              ybar=ybar,
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
                      lambda=lambda
                     )
  A = Prod$A
  a = Prod$a
  
  # Amenities
  AM = living_amenities(theta=theta,
                 N=N,
                 L_i=L_i,
                 W_i=W_i,
                 Q=Q,
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
  U = (sumDims(array_operator(u, array_operator(U_i, theta, '^'), '*'),1))^(1/theta)

  return(list(A=A, a=a, u=u, B=B, b=b, w=w, varphi=varphi, U=U, Q_norm=Q_norm))
}


#' Function to solve counterfactuals.
#'
#' @param N Integer - Number of locations.
#' @param L_i Nx1 array - Number of residents in each location
#' @param L_j Nx1 array - Number of workers in each location
#' @param t_ij NxN matrix - Travel times across locations
#' @param varphi Nx1 array - Density of development
#' @param K Nx1 array - Land supply
#' @param a Nx1 array - Total Factor Productivity in each location
#' @param b Nx1 array - Vector of amenities in each location
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' @param endo_Lr 0: workers don't reallocate 1: workers reallocate
#' @param alpha Float - Exp. share in consumption, 1-alpha exp. share in housing
#' @param beta Float - Output elasticity with respect to labor
#' @param theta Float - Commuting and migration elasticity.
#' @param mu Float - Floorspace prod function: output elasticity wrt capital
#' @param delta Float - Decay parameter agglomeration force
#' @param lambda Float - agglomeration externality
#' @param rho Float - decay parameter for amenities
#' @param eta Float - amenity externality
#' @param epsilon Float - Parameter that transforms travel times to commuting costs
#' @param w_eq Nx1 array - Initial vector of wages
#' @param u_eq Nx1 array - Initial vector of welfare
#' @param Q_eq Nx1 array - Initial price for floorspace
#' @param ttheta_eq Nx1 array - Share of floorspace used commercially 
#' @param zeta Float - convergence parameter
#' 
#' @return Counterfactual values.
#' @export
#'
#' @examples
solveModel = function(N,
                      L_i,
                      L_j,
                      varphi,
                      t_ij,
                      K,
                      a,
                      b,
                      maxiter,
                      endo_Lr,
                      alpha,
                      beta,
                      theta,
                      mu,
                      delta,
                      lambda,
                      rho,
                      eta,
                      epsilon,
                      w_eq,
                      u_eq,
                      Q_eq,
                      theta_eq
                      ){

  # Formatting of input data
  D = commuting_matrix(t_ij=t_ij, epsilon = epsilon)
  tau = D$tau
  L_i = array(unlist(L_i),dim(L_i))
  L_j = array(unlist(L_j),dim(L_j))
  K = array(unlist(K), dim(K))

  # Settings
  outerdiff = Inf;
  w = w_eq;
  u = u_eq;
  Q = Q_eq;
  ttheta = ttheta_eq
  iter = 0;
  zeta_init = zeta;

  while(outerdiff>tol & iter < maxiter){
    # 1) Labor supply equation
    w_tr = aperm(array(w, dim=c(N,1)), c(2,1));
    rep_w_tr = kronecker(w_tr^theta, array(1, dim=c(N, 1)));
    # Constructing employment shares
    w_tr_tau = array_operator(w_tr^theta, tau^(-theta), '*');
    lambda_ij_i = array_operator(w_tr_tau, sumDims2(w_tr_tau,2), '/');
    W_i = (sumDims2(w_tr_tau,2))^(1/theta);
    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    L_ij = array_operator(L_i, lambda_ij_i, '*')
    L_j = sumDims2(L_ij, 1)
    L = sum(L_i)
    
    # 2 average income
    av_income = av_income(lambda_ij_i=lambda_ij_i,w_tr = w_tr)
    ybar = av_income$y_bar
    
    # 3 Total floorspace
    FS = array_operator(varphi,K^(1-mu),"*")
    
    # 4 Agglomeration externalities
    L_j_dens = (array_operator(L_j, K, '/'));
    L_j_dens_per = aperm(array(L_j_dens, dim=c(N,1)), c(2,1));
    L_j_dens_rep = kronecker(L_j_dens_per, array(1, dim=c(N, 1)));
    Upsilon = sumDims2(array_operator(exp(-delta*t_ij), L_j_dens_rep, '*'), 2);    
    A = array_operator(a, Upsilon, '*')
    
    # 5 Amenities
    L_i_dens = (array_operator(L_i, K, '/'));
    L_i_dens_per = aperm(array(L_i_dens, dim=c(N,1)), c(2,1));
    L_i_dens_rep = kronecker(L_i_dens_per, array(1, dim=c(N, 1)));
    Omega = sumDims2(array_operator(exp(-rho*t_ij), L_i_dens_rep, '*'), 2);
    B = array_operator(b, Omega,'*')
    
    # 6 Residents, probabilities, and welfare
    u =  array_operator(array_operator(W_i, Q^(1-alpha), '/'), B, '*')
    U = sum(u^theta)
    lambda_i = (u^theta)/U
    L_i_upd = L*lambda_i
    
    # 7 Total output by location
    FS_f = array_operator(ttheta,array_operator(varphi, K^(1-mu), '*'), '*')
    Y = array_operator(A, array_operator(L_j^beta, FS_f^beta, '*'), '*')
    Q_upd1 = (1-beta)*array_operator(Y,FS_f)
    w_upd = beta*array_operator(Y, L_j, '/')

    # 8 Housing prices
    FS_r = array_operator((1-ttheta), array_operator(varphi, K^(1-mu), '*'), '*')
    X = array_operator(ybar, L_i, '*')
    Q_upd2 = (1-alpha)*array_operator(X, FS_r, '/')
    Q_upd = Q_upd1*(a>0) + Q_upd2*(a==0 & b>0)
    
    # 9 Share of commercial floorspace
    LP = array_operator(Q_upd1, array_operator(varphi, K^(1-mu), '*'), '*')
    ttheta_upd = (1-beta)*array_operator(Y, LP, '/')
    ttheta_upd = (b==0)+ttheta_upd*(b>0)
    
    # 10 Calculating the main differences
    z_w = array_operator(w, w_upd, '-')
    z_L = array_operator(L_i, L_i_upd, '-')
    z_Q = array_operator(Q, Q_upd, '-')
    z_theta = array_operator(ttheta, ttheta_upd, '-')
    outerdiff = max(max(abs(z_w)), max(abs(z_L)), max(abs(z_Q)), max(abs(z_theta)))
    iter = iter+1
    
    # 11 New vector of variables
    L_i = zeta*L_i + (1-zeta)*L_i_upd
    Q = zeta*Q + (1-zeta)*Q_upd
    w = zeta*w + (1-zeta)*w_upd
    ttheta = zeta*ttheta + (1-zeta)*ttheta_upd
    
  }
  
  return(list(w=w, W_i=W_i, B=B, A=A, Q=Q, lambda_ij_i=lambda_ij_i, L_i=L_i, L_j=L_j,
              ybar=ybar, ttheta=ttheta, u=u, U=U))
}

