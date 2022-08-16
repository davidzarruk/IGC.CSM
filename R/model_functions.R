
#' Function to compute equilibrium wages that make the model labor in every
#' location in every sector equal to the observed data. It finds the w's
#' such that equation (3.2) holds.
#'
#' @param N Integer - Number of locations.
#' @param S Integer - Number of sectors in the economy.
#' @param w_init Wage initialization.
#' @param theta1 Float - Commuting elasticity.
#' @param tau NxN matrix - Commuting costs between all possible locations.
#' @param lambda_is_i NxS matrix - Proportion of residents in each location in
#'     each sector. Rows add up to 1.
#' @param lambda_i Nx1 matrix - Number of residents in each location.
#' @param L_j_data NxS matrix - Total amount of workers in each location in each
#'     sector.
#' @param nu_init Float - Convergence parameter to update wages.
#'     Default nu=0.005.
#' @param tol Float - Maximum tolerable error for estimating total labor.
#'     Default tol=10^-6.
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' @param nu_mult Vector - Multipliers of nu.
#'     Default is empty.
#' @param nu_intervals Vector - Intervals to use multipliers.
#'     Default is empty.
#'
#' @return A list with equilibrium wages and probability of workers in each
#'     location working in every other location.
#' @export
#'
#' @examples
labor_productivity = function(N,
                              S,
                              w_init,
                              theta1,
                              tau,
                              lambda_is_i,
                              lambda_i,
                              L_j_data,
                              nu_init=0.005,
                              tol=10^-6,
                              maxiter=1000,
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
    w_tr = aperm(array(w, dim=c(N,S,1)), c(3,1,2))
    reptau = kronecker(tau, array(1, dim=c(1, 1, S)))

    # Constructing probabilities
    w_tr_reptau = array_operator(w_tr^theta1, reptau^(-theta1), '*')
    lambda_ijs_is = array_operator(w_tr_reptau, sumDims(w_tr_reptau, 2), '/')
    lambda_is_i_p = aperm(kronecker(lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));

    # Labor is equal to probabilities * total number of residents * proportion of workers in each sector.
    l_prod = array_operator(lambda_i, array_operator(lambda_is_i_p, lambda_ijs_is, '*'), '*')
    L_j_p = sumDims(l_prod, 1)
    L_j = aperm(L_j_p, c(2, 3, 1));

    z_L = array_operator(L_j_data, L_j, '-')
    w = array_operator(w, 1+nu*(z_L/L_j), '*');
    w = w/w[1,1,1];
    outerdiff = max(abs(z_L))

    iter = iter+1;

    if(iter>100){
      nu = nu_init*10;
    } else if(iter>1000){
      nu = nu_init*50;
    } else if(iter>10000){
      nu = nu_init*100;
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

  return(list(w=w, w_tr=w_tr, lambda_ijs_is=lambda_ijs_is, l_prod=l_prod))
}


#' Function to transform travel times into iceberg commuting costs
#' @param t_ij NxN matrix - Travel time matrix across locations
#' @param delta Float - Parameter that transforms travel times to commuting costs
#' 
#' @return A NxN matrix of commuting costs
#' 
#' @export
#' 
#'  @examples 
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


#' Function to recover the amenity distribution of working in each location.
#' Finds B's that make the probability of working in each sector in equation
#' (3.2) equal to the probabilities observed in the data.
#'
#' @param theta1 Float - Commuting elasticity.
#' @param kappa1 Float - Parameter that governs the reallocation of workers
#'     across sectors. It measures how easy it is for workers to substitute
#'     jobs across sectors. In the case in which it tends to infity, the model
#'     replicates the case in which there are no frictions in the labor market,
#'     and workers are perfectly mobile across sectors. The case in which it
#'     tends to one replicates the specific factor model in which workers are
#'     only productive in one sector. Then, they do not switch jobs across
#'     industries.
#' @param lambda_is_i NxS matrix - Proportion of residents in each location in
#'     each sector. Rows add up to 1.
#' @param B_init Matrix - Initialization parameter for the amenity distribution
#'     of work.
#' @param w_tr Matrix - Equilibrium wages
#' @param reptau NxN matrix - Commuting costs between all possible locations.
#' @param tol Float - Maximum tolerable error for estimating total labor.
#'     Default tol=10^-6.
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#'
#' @return
#' @export
#'
#' @examples
amenities = function(theta1,
                     kappa1,
                     lambda_is_i,
                     B_init,
                     w_tr,
                     reptau,
                     tol,
                     maxiter){

  B = B_init;
  outerdiff = Inf;
  iter = 0;

  while(outerdiff>tol & iter<maxiter){
    # Amenities
    W_is_aux_1 = sumDims(array_operator(w_tr^theta1, reptau^(-theta1), '*'), 2)
    W_is_aux_1 = W_is_aux_1^(1/theta1)
    W_is = aperm(W_is_aux_1, c(1,3,2))

    B_W_is = array_operator(B, W_is, '*')^(kappa1)

    lambda_is_i_m = array_operator(B_W_is, sumDims(B_W_is, 2), '/');
    lambda_is_i_m_B = array_operator(lambda_is_i_m, B, '/');

    B_upd = array_operator(lambda_is_i, lambda_is_i_m_B, '/');
    B_upd = array_operator(B_upd, B_upd[,1,], '/');

    z_B = array_operator(lambda_is_i, lambda_is_i_m, '-');
    diff_B = max(abs(z_B));
    B = 0.2*B_upd+0.8*kronecker(B, array(1, dim=c(1,1,1)));
    B = array_operator(B, B[,1,1], '/');
    outerdiff = diff_B

    iter = iter+1;

    print(outerdiff)
  }

  return(list(W_is=W_is, B=B))
}

#' Computes average income in each location, which is the weighted average of
#' the income of the people living in the location.
#'
#' @param lambda_is_i NxS matrix - Proportion of residents in each location in
#'     each sector. Rows add up to 1.
#' @param B Amenity parameters of each location.
#' @param W_is
#' @param kappa1 Float - Parameter that governs the reallocation of workers
#'     across sectors. It measures how easy it is for workers to substitute
#'     jobs across sectors. In the case in which it tends to infinity, the model
#'     replicates the case in which there are no frictions in the labor market,
#'     and workers are perfectly mobile across sectors. The case in which it
#'     tends to one replicates the specific factor model in which workers are
#'     only productive in one sector. Then, they do not switch jobs across
#'     industries.
#' @param lambda_ijs_is NxNxS matrix - Probability of individuals in each
#'     location of working in each location and each sector.
#' @param L_bar Float - Total number of people in the city.
#' @param w NxS - Wages in each location in each sector.
#' @param lambda_i Nx1 matrix - Number of residents in each location.
#'
#' @return
#' @export
#'
#' @examples
av_income = function(lambda_is_i,
                     B,
                     W_is,
                     kappa1,
                     lambda_ijs_is,
                     L_bar,
                     w,
                     lambda_i){

  lambda_is_i_p = aperm(kronecker(lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));
  y_bar = (sumDims(array_operator(B, W_is, '*')^(kappa1), 2))^(1/kappa1);
  L_j_eff_p = sumDims(array_operator(y_bar, array_operator(lambda_i, array_operator(lambda_is_i_p, array_operator(lambda_ijs_is, L_bar, '*'), '*'), '*'), '*'),1);
  L_j_eff_w = aperm(L_j_eff_p, c(2, 3, 1));
  L_j_eff = array_operator(L_j_eff_w, w, '/');

  return(list(lambda_is_i_p=lambda_is_i_p,
              L_j_eff_p=L_j_eff_p,
              L_j_eff_w=L_j_eff_w,
              L_j_eff=L_j_eff,
              y_bar=y_bar))
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
av_income_simple = function(lambda_ij_i,
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



#' Computes residential and commercial floorspace supply and equilibrium prices.
#'
#' @param alpha1 Float - Utility parameter that determines preferences for
#'     consumption.
#' @param H_bar_rest Float -
#' @param H_bar Nx1 array - Land supply.
#' @param w NxS - Wages in each location in each sector.
#' @param L_j_eff
#' @param y_bar - Average income in each location.
#' @param lambda_i Nx1 matrix - Number of residents in each location.
#' @param L_bar Float - Total number of people in the city.
#' @param beta1 Float - Cobb-Douglas parameter for the production function on
#'     each sector.
#'
#' @return
#' @export
#'
#' @examples
floorspace_supply = function(alpha1,
                             H_bar_rest,
                             H_bar,
                             w,
                             L_j_eff,
                             y_bar,
                             lambda_i,
                             L_bar,
                             beta1){

  H_r = pmin(H_bar_rest, H_bar);
  H_f = pmin(H_bar_rest, H_bar);
  beta_tilde_sector = aperm(array_operator((1-beta1), beta1, '/'), c(1, 3, 2));
  q = (array_operator(sumDims(array_operator(beta_tilde_sector, array_operator(w, L_j_eff, '*'), '*'),2), H_f, '/'));
  r = ((1-alpha1)/alpha1)*(array_operator(array_operator(y_bar, lambda_i*L_bar, '*')+array_operator(q, H_f, '*'), H_r, '/'));

  return(list(beta_tilde_sector=beta_tilde_sector, q=q, r=r, H_r=H_r, H_f=H_f))
}

#' Computes number of firms
#'
#' @param N Int - Number of locations.
#' @param beta1 Float - Cobb-Douglas parameter for the production function on
#'     each sector.
#' @param beta_tilde_sector Float - (1-beta)/beta
#' @param w NxS - Wages in each location in each sector.
#' @param L_j_eff
#' @param sigma1 - elasticity of subsitution across varieties within each sector.
#' @param q - Commercial floorspace prices.
#'
#' @return
#' @export
#'
#' @examples
number_firms = function(N,
                        beta1,
                        beta_tilde_sector,
                        w,
                        L_j_eff,
                        sigma1,
                        q,
                        F){

  beta_tilde = aperm(array_operator(array_operator(beta1, (-beta1), '^'), (array_operator((1-beta1), (-(1-beta1)), '^')), '*'), c(1, 3, 2));
  H_f_sector = array_operator((array_operator(kronecker(beta_tilde_sector, array(1, dim=c(N,1,1))), sumDims(array_operator(w, L_j_eff, '*'), 2), '*')), q, '/');
  sigma_cons = aperm(sigma1, c(1, 3, 2));
  beta_factor = aperm(beta1, c(1, 3, 2));
  M = array_operator(array_operator(array_operator(beta_tilde, array_operator(L_j_eff, (beta_factor), '^'), '*'), array_operator(H_f_sector,(1-beta_factor), '^'), '*'), array_operator(sigma_cons, F, '*'), '/');

  return(list(beta_factor=beta_factor, M=M, sigma_cons=sigma_cons, beta_tilde=beta_tilde))
}

#' Computes prices charged by firms in each location to consumers in each
#'     location.
#'
#' @param N Int - Number of locations.
#' @param sigma1 - elasticity of subsitution across varieties within each sector.
#' @param beta_factor
#' @param w NxS - Wages in each location in each sector.
#' @param q - Commercial floorspace prices.
#' @param A
#'
#' @return
#' @export
#'
#' @examples
price_location = function(N,
                          sigma1,
                          beta_factor,
                          w,
                          q,
                          A){

  markup = aperm(array_operator(sigma1, (sigma1-1), '/'), c(1, 3, 2));
  inv_sigma = aperm((1/(1-sigma1)), c(1, 3, 2));
  p_j = array_operator(array_operator(array_operator(markup, array_operator(w, beta_factor, '^'), '*'), array_operator(q, kronecker(1-beta_factor, array(1, dim=c(N,1,1))), '^'), '*'), A, '/');

  return(list(p_j=p_j, inv_sigma=inv_sigma, markup=markup))
}

#' Title
#'
#' @param sigma1
#' @param xi1
#' @param M
#' @param p_j
#' @param sigma1_tr
#' @param inv_sigma
#'
#' @return
#' @export
#'
#' @examples
agg_price_index = function(sigma1,
                           xi1,
                           M,
                           p_j,
                           sigma1_tr,
                           inv_sigma){

  sigma1_tr = aperm(sigma1, c(1, 3, 2));
  P_s = array_operator((sumDims(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'),1)), (inv_sigma), '^');
  pi_js = array_operator(array_operator(M, array_operator(p_j, (1-sigma1_tr), '^'), '*'), (array_operator(P_s, (1-sigma1_tr), '^')), '/');
  P = array_operator(sumDims(array_operator(P_s, (1-xi1), '^'),2), (1/(1-xi1)), '^');
  pi_s = array_operator(array_operator(P_s, (1-xi1), '^'), (array_operator(P, (1-xi1), '^')), '/');

  return(list(pi_js=pi_js, pi_s=pi_s, P=P, P_s=P_s, sigma1_tr=sigma1_tr))
}

total_sales = function(N,
                       y_bar,
                       lambda_i,
                       L_bar,
                       r,
                       H_r,
                       q,
                       H_f,
                       alpha1,
                       pi_js,
                       pi_s){

  X = sumDims(array_operator(y_bar, array_operator(lambda_i, kronecker(L_bar, array(1, dim=c(N))), '*'), '*') + array_operator(r, H_r, '*') +array_operator(q, H_f, '*'), 1);
  Y_js = array_operator(alpha1, array_operator(pi_js, array_operator(pi_s, X, '*'), '*'), '*');

  return(list(X=X, Y_js=Y_js))
}


#' Title
#'
#' @param sigma_cons
#' @param w
#' @param L_j_eff
#' @param beta_factor
#' @param Y_js
#' @param A
#' @param inv_sigma
#'
#' @return
#' @export
#'
#' @examples
new_labor_demand = function(sigma_cons,
                            w,
                            L_j_eff,
                            beta_factor,
                            Y_js,
                            A,
                            inv_sigma){

  inv_sigma = array_operator(sigma_cons-1, -1, '^');
  LS = array_operator(w, L_j_eff, '*');
  LD = array_operator(beta_factor, Y_js, '*');
  LD_A = array_operator(LD, array_operator(A, (sigma_cons-1), '^'), '/');
  A_prime = array_operator(array_operator(LS, LD_A, '/'), inv_sigma, '^');

  return(list(inv_sigma=inv_sigma, LS=LS, LD=LD, LD_A=LD_A, A_prime=A_prime))
}


#' Title
#'
#' @param alpha1
#' @param beta1
#' @param eta1
#' @param kappa1
#' @param sigma1
#' @param xi1
#' @param L_bar
#' @param H_bar
#' @param H_bar_rest
#' @param N
#' @param S
#' @param A_init
#' @param lambda_i
#' @param lambda_is_i
#' @param w
#' @param W_is
#' @param B
#' @param lambda_ijs_is
#' @param zeta_init
#' @param zeta
#' @param tol
#' @param maxiter
#' @param zeta_mult
#' @param zeta_intervals
#'
#' @return
#' @export
#'
#' @examples
eq_quantities = function(alpha1,
                         beta1,
                         eta1,
                         kappa1,
                         sigma1,
                         xi1,
                         L_bar,
                         H_bar,
                         H_bar_rest,
                         N,
                         S,
                         A_init,
                         F,
                         lambda_i,
                         lambda_is_i,
                         w,
                         W_is,
                         B,
                         lambda_ijs_is,
                         zeta_init = zeta_init,
                         zeta,
                         tol,
                         maxiter,
                         zeta_mult = c(),
                         zeta_intervals = c()){
  outerdiff = Inf;
  A = A_init;
  iter = 0;

  # 2) Average income in each location
  av_inc = av_income(lambda_is_i, B, W_is, kappa1, lambda_ijs_is, L_bar, w, lambda_i)
  L_j_eff = av_inc$L_j_eff
  y_bar = av_inc$y_bar

  # 3) Floorspace supply
  floorspace_sup = floorspace_supply(alpha1, H_bar_rest, H_bar, w, L_j_eff, y_bar, lambda_i, L_bar, beta1)
  beta_tilde_sector = floorspace_sup$beta_tilde_sector
  q = floorspace_sup$q
  r = floorspace_sup$r
  H_r = floorspace_sup$H_r
  H_f = floorspace_sup$H_f

  # 5) Number of firms by location
  number_fir = number_firms(N, beta1, beta_tilde_sector, w, L_j_eff, sigma1, q, F)
  beta_factor = number_fir$beta_factor
  M = number_fir$M
  sigma_cons = number_fir$sigma_cons

  while(outerdiff>tol & iter < maxiter){
    # 6) Price in each location
    price_loc = price_location(N, sigma1, beta_factor, w, q, A)
    p_j = price_loc$p_j
    inv_sigma = price_loc$inv_sigma

    # 7) Aggregate price index and expenditure share in the variety of each location
    agg_price_ind = agg_price_index(sigma1, xi1, M, p_j, sigma1_tr, inv_sigma)
    pi_js = agg_price_ind$pi_js
    pi_s = agg_price_ind$pi_s
    P = agg_price_ind$P

    # 8) Total sales in each location
    total_sal = total_sales(N, y_bar, lambda_i, L_bar, r, H_r, q, H_f, alpha1, pi_js, pi_s)
    Y_js = total_sal$Y_js

    # 9) New labor demand
    new_labor_dem = new_labor_demand(sigma_cons, w, L_j_eff, beta_factor, Y_js, A, inv_sigma)
    A_prime = new_labor_dem$A_prime

    # 10) Update
    z_A = array_operator(A, A_prime, '-');
    A = array_operator(zeta*A_prime, (1-zeta)*A, '+');
    outerdiff = max(abs(z_A))
    iter = iter + 1;

    print(outerdiff)

    if(outerdiff<500){
      zeta = zeta_init*10;
    }
    if(outerdiff<300){
      zeta = zeta_init*50;
    }
    if(outerdiff<100){
      zeta = zeta_init*100;
    }
    if(outerdiff<20){
      zeta = zeta_init*500;
    }
    if(outerdiff<10){
      zeta = zeta_init*500;
    }
    if(outerdiff<1){
      zeta = zeta_init*500;
    }
    if(iter>100000){
      break
    }

    # if(length(zeta_mult) > 0){
    #   if(length(zeta_mult) != length(zeta_intervals)){
    #     stop(
    #       "The array zeta_mult must have the same size as the array zeta_intervals."
    #     )
    #   } else{
    #     iiter = 1
    #     while(iter > zeta_intervals[iiter] & iiter <= length(zeta_intervals)){
    #       iiter = iiter+1
    #     }
    #     zeta=zeta_init*zeta_mult[iiter]
    #   }
    # }
  }

  # Amenity parameters
  num_U_i = array_operator(sumDims(array_operator(W_is, kappa1, '^'),2), (1/kappa1), '^');
  den_U_i = array_operator(array_operator(r, (1-alpha1), '^'), array_operator(P, alpha1, '^'), '*');
  U_i = array_operator(num_U_i, den_U_i, '/');

  return(list(U_i=U_i, A=A))
}

#' Function to estimate amenity parameters of locations where users live.
#'
#' @param eta1 Float - Parameter that governs the reallocation of workers across
#'     locations in the city. This parameter measures how sensible are migration
#'     flows within the city to changes in real income.
#' @param lambda_i Nx1 matrix - Probability of living in each location.
#' @param u_init Nx1 matrix - initialization matrix for amenity parameters.
#' @param U_i Matrix - Real income in each location.
#' @param zeta Float - Parameter to update amenity parameters in each iteration.
#' @param tol Float - Maximum tolerable error for estimating total labor.
#'     Default tol=10^-6.
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#'
#' @return Matrix with the amenity distribution of living in each location.
#' @export
#'
#' @examples
living_amenities = function(eta1,
                              lambda_i,
                              u_init,
                              U_i,
                              zeta,
                              tol,
                              maxiter){
  outerdiff = Inf;
  u = u_init;
  iter = 0;

  while(outerdiff>tol & iter < maxiter){
    lambda_i_model = array_operator(array_operator(u, array_operator(U_i, eta1, '^'), '*'), (sumDims(array_operator(u, array_operator(U_i, eta1, '^'), '*'),1)), '/');
    lambda_i_m_u = array_operator(lambda_i_model, u, '/');
    u_prime = array_operator(lambda_i, lambda_i_m_u, '/');

    u = array_operator(zeta*u_prime, (1-zeta)*u, '+');
    u = array_operator(u, u[1,1,1], '/');
    outerdiff = max(abs(array_operator(lambda_i_model, lambda_i, '-')))
    print(outerdiff)
  }
  return(list(u=u))
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
living_amenities_simple = function(theta,
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

#' Function to invert model, so amenities and wages are chosen to match model
#' to data.
#'
#' @param N Integer - Number of locations.
#' @param S Integer - Number of sectors in the economy.
#' @param L_bar TODO
#' @param H_bar TODO
#' @param H_bar_rest TODO
#' @param tau NxN matrix - Commuting costs between all possible locations.
#' @param lambda_i Nx1 matrix - Number of residents in each location.
#' @param lambda_is_i NxS matrix - Proportion of residents in each location in
#'     each sector. Rows add up to 1.
#' @param L_j_data NxS matrix - Total amount of workers in each location in each
#'     sector.
#' @param zeta Int - Convergence parameter for amenities
#' @param z_init TODO
#' @param tol TODO
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' @param alpha1 Float - Utility parameter that determines preferences for
#'     consumption.
#' @param beta0 TODO
#' @param theta1 Float - Commuting elasticity.
#' @param eta1 TODO
#' @param kappa1 Float - Parameter that governs the reallocation of workers
#'     across sectors. It measures how easy it is for workers to substitute
#'     jobs across sectors. In the case in which it tends to infity, the model
#'     replicates the case in which there are no frictions in the labor market,
#'     and workers are perfectly mobile across sectors. The case in which it
#'     tends to one replicates the specific factor model in which workers are
#'     only productive in one sector. Then, they do not switch jobs across
#'     industries.
#' @param sigma0 TODO
#' @param xi1 TODO
#' @param nu_init TODO
#' @param nu_intervals TODO
#' @param nu_mult TODO
#' @param zeta_intervals TODO
#' @param zeta_mult TODO
#'
#' @return Equilibrium values.
#' @export
#'
#' @examples
inversionModel_Eff = function(N,
                              S,
                              L_bar,
                              H_bar,
                              H_bar_rest,
                              tau,
                              lambda_i,
                              lambda_is_i,
                              L_j_data,
                              zeta,
                              z_init,
                              tol=10^-6,
                              maxiter=1000,
                              alpha1=0.7,
                              beta0=0.7,
                              theta1=7,
                              eta1=1.5,
                              kappa1=2,
                              sigma0=6,
                              xi1=1.8,
                              nu_init=0.005,
                              F=1,
                              nu_intervals = c(),
                              nu_mult = c(),
                              zeta_intervals = c(),
                              zeta_mult = c()){

  # Formatting of input data
  tau = array(unlist(tau), dim(tau))
  H_bar = array(unlist(H_bar), dim(H_bar))
  L_j_data = array(unlist(L_j_data), dim(L_j_data))
  lambda_is_i = array(unlist(lambda_is_i), dim(lambda_is_i))
  lambda_i = array(unlist(lambda_i), dim(lambda_i))

  # Parameters
  beta1 = beta0*array(1, dim=c(1,1,S))
  sigma1 = sigma0*array(1, dim=c(1,1,S))

  # Initialization
  A_init=array(1, dim=c(N,S))
  u_init=array(1, dim=c(N,1))
  B_init=array(1, dim=c(N,S))
  w_init=array(1, dim=c(N,S))

  # Labor productivity
  LP = labor_productivity(N=N,
                          S=S,
                          w_init=w_init,
                          theta1=theta1,
                          tau=tau,
                          lambda_is_i=lambda_is_i,
                          lambda_i=lambda_i,
                          L_j_data=L_j_data,
                          nu_init=nu_init,
                          tol=tol,
                          maxiter=maxiter,
                          nu_intervals = nu_intervals,
                          nu_mult = nu_mult)

  # Equilibrium wages
  w = LP$w
  w_tr = LP$w_tr

  reptau = kronecker(tau, array(1, dim=c(1, 1, S)))

  # Amenities
  AM = amenities(theta1=theta1,
                 kappa1=kappa1,
                 lambda_is_i=lambda_is_i,
                 B_init=B_init,
                 w_tr,
                 reptau,
                 tol,
                 maxiter=maxiter)

  W_is = AM$W_is
  B = AM$B

  # Other equilibrium quantities
  print("eq quantities starting")
  EQ = eq_quantities(alpha1=alpha1,
                     beta1=beta1,
                     eta1=eta1,
                     kappa1=kappa1,
                     sigma1=sigma1,
                     xi1=xi1,
                     L_bar=L_bar,
                     H_bar=H_bar,
                     H_bar_rest=H_bar_rest,
                     N=N,
                     S=S,
                     A_init=A_init,
                     F=F,
                     lambda_i=lambda_i,
                     lambda_is_i=lambda_is_i,
                     w=w,
                     W_is=W_is,
                     B=B,
                     lambda_ijs_is = LP$lambda_ijs_is,
                     zeta_init = z_init,
                     zeta=zeta,
                     tol=tol,
                     maxiter=maxiter,
                     zeta_mult = zeta_mult,
                     zeta_intervals = zeta_intervals)

  U_i = EQ$U_i
  A = EQ$A

  LEQ = living_amenities(eta1=eta1,
                         lambda_i=lambda_i,
                         u_init=u_init,
                         U_i=U_i,
                         zeta=zeta,
                         tol=tol,
                         maxiter=maxiter)
  u = LEQ$u

  # Save and export
  U = (sumDims(array_operator(u, array_operator(U_i, eta1, '^'), '*'),1))^(1/eta1)

  return(list(A=A, u=u, B=B, w=w))
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
  Inc = av_income_simple(lambda_ij_i=lambda_ij_i,
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
  AM = living_amenities_simple(theta=theta,
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
#' @param S Integer - Number of sectors in the economy.
#' @param L_bar TODO
#' @param H_bar TODO
#' @param H_bar_rest TODO
#' @param tau NxN matrix - Commuting costs between all possible locations.
#' @param A TODO
#' @param u_eq TODO
#' @param B TODO
#' @param w_eq TODO
#' @param lambda_i Nx1 matrix - Number of residents in each location.
#' @param lambda_is_i NxS matrix - Proportion of residents in each location in
#'     each sector. Rows add up to 1.
#' @param zeta TODO
#' @param tol TODO
#' @param maxiter Integer - Maximum number of iterations for convergence.
#'     Default maxiter=1000.
#' @param endo_Lr TODO
#' @param alpha1 Float - Utility parameter that determines preferences for
#'     consumption.
#' @param beta0 TODO
#' @param theta1 Float - Commuting elasticity.
#' @param eta1 TODO
#' @param kappa1 Float - Parameter that governs the reallocation of workers
#'     across sectors. It measures how easy it is for workers to substitute
#'     jobs across sectors. In the case in which it tends to infity, the model
#'     replicates the case in which there are no frictions in the labor market,
#'     and workers are perfectly mobile across sectors. The case in which it
#'     tends to one replicates the specific factor model in which workers are
#'     only productive in one sector. Then, they do not switch jobs across
#'     industries.
#' @param xi1 TODO
#' @param sigma0 TODO
#'
#' @return Counterfactual values.
#' @export
#'
#' @examples
solveModel1_Eff = function(N,
                           S,
                           L_bar,
                           H_bar,
                           H_bar_rest,
                           tau,
                           A,
                           u_eq,
                           B,
                           w_eq,
                           lambda_i,
                           lambda_is_i,
                           zeta,
                           tol=10^-6,
                           maxiter=1000,
                           endo_Lr=1,
                           alpha1=0.7,
                           beta0=0.7,
                           theta1=7,
                           eta1=1.5,
                           kappa1=2,
                           xi1=1.8,
                           sigma0=6,
                           F=1){

  # Formatting of input data
  tau = array(unlist(tau), dim(tau))
  H_bar = array(unlist(H_bar), dim(H_bar))
  lambda_is_i = array(unlist(lambda_is_i), dim(lambda_is_i))
  lambda_i = array(unlist(lambda_i), dim(lambda_i))

  # Parameters
  beta1 = beta0*array(1, dim=c(1,1,S))
  sigma1 = sigma0*array(1, dim=c(1,1,S))

  # Settings
  outerdiff = Inf;
  w = w_eq;
  u = u_eq;
  iter = 0;
  zeta_init = zeta;

  while(outerdiff>tol & iter < maxiter){
    # 1) Labor supply - Indirect utility
    w_tr = aperm(array(w, dim=c(N,S,1)), c(3,1,2))
    reptau = kronecker(tau, array(1, dim=c(1, 1, S)))

    w_tr_reptau = array_operator(w_tr^theta1, reptau^(-theta1), '*')
    lambda_ijs_is = array_operator(w_tr_reptau, sumDims(w_tr_reptau, 2), '/')

    W_is_aux_1 = sumDims(w_tr_reptau,2);
    W_is_aux = W_is_aux_1^(1/theta1);
    W_is = aperm(W_is_aux, c(1,3,2));

    B_W_is = array_operator(B, W_is, '*')^(kappa1)
    lambda_is_i = array_operator(B_W_is, sumDims(B_W_is,2), '/');

    lambda_is_i_p = aperm(kronecker(lambda_is_i, array(1, dim=c(1,1,1))), c(1, 3, 2));

    l_prod = array_operator(lambda_i, array_operator(lambda_is_i_p, lambda_ijs_is, '*'), '*')
    L_j_p = sumDims(l_prod, 1)
    L_j = aperm(L_j_p, c(2, 3, 1));

    # 2) Average income in each location
    av_inc = av_income(lambda_is_i=lambda_is_i,
                       B=B,
                       W_is=W_is,
                       kappa1=kappa1,
                       lambda_ijs_is=lambda_ijs_is,
                       L_bar=L_bar,
                       w=w,
                       lambda_i=lambda_i);
    L_j_eff=av_inc$L_j_eff
    y_bar=av_inc$y_bar

    # 3) Floorspace supply
    floorspace = floorspace_supply(alpha1, H_bar_rest, H_bar, w, L_j_eff, y_bar, lambda_i, L_bar, beta1)
    H_r = floorspace$H_r
    H_f = floorspace$H_f
    beta_tilde_sector = floorspace$beta_tilde_sector
    q = floorspace$q
    r = floorspace$r

    # 5) Number of firms by location
    num_firms = number_firms(N, beta1, beta_tilde_sector, w, L_j_eff, sigma1, q, F)
    beta_factor=num_firms$beta_factor
    M=num_firms$M
    sigma_cons=num_firms$sigma_cons

    H_f_sector = array_operator(array_operator(kronecker(beta_tilde_sector, array(1, dim=c(N,1,1))), sumDims(array_operator(w, L_j_eff, '*'),2), '*'), q, '/');

    Q_js = array_operator(array_operator(L_j_eff, beta_factor, '^'), array_operator(H_f_sector, 1-beta_factor, '^'), '*');
    elast_demand = aperm(array_operator(sigma1, sigma1-1, '/'), c(1,3,2));
    inv_elast_demand = aperm(array_operator(sigma1-1, sigma1, '/'), c(1,3,2));

    Q_s = array_operator(sumDims(array_operator(Q_js, elast_demand, '^'), 1), inv_elast_demand, '^');
    Q = array_operator(sumDims(array_operator(Q_s, (xi1-1)/xi1, '^'),2), xi1/(xi1-1), '^');

    # 6) Price in each location
    price_loc = price_location(N, sigma1, beta_factor, w, q, A)
    p_j = price_loc$p_j
    inv_sigma = price_loc$inv_sigma

    # 7) Aggregate price index and expenditure share in the variety of each location
    agg_price_ind = agg_price_index(sigma1, xi1, M, p_j, sigma1_tr, inv_sigma)
    pi_js = agg_price_ind$pi_js
    pi_s = agg_price_ind$pi_s
    P = agg_price_ind$P

    # 8) Total sales in each location
    total_sal = total_sales(N, y_bar, lambda_i, L_bar, r, H_r, q, H_f, alpha1, pi_js, pi_s)
    Y_js = total_sal$Y_js

    # 9) New labor demand and wages
    new_labor_dem = new_labor_demand(sigma_cons, w, L_j_eff, beta_factor, Y_js, A, inv_sigma)
    LD = new_labor_dem$LD
    w_prime = array_operator(LD, L_j_eff, '/');
    w_prime = array_operator(w_prime, w_prime[1,1,1], '/');

    # 10) Update
    w = zeta*w_prime + (1-zeta)*w;
    w = array_operator(w, w[1,1,1], '/');

    w_diff = max(abs(w - w_prime));

    num_U_i = (sumDims(array_operator(W_is, kappa1, '^'),2))^(1/kappa1);
    den_U_i = array_operator((r^(1-alpha1)), (P^(alpha1)), '*');
    U_i = array_operator(num_U_i, den_U_i, '/');
    U = array_operator(sumDims(array_operator(u, U_i^eta1, '*'),1), 1/eta1, '^');

    if(endo_Lr==1){
      lambda_i_upd = array_operator(array_operator(u, U_i^eta1, '*'), sumDims(array_operator(u, U_i^eta1, '*'),1), '/');
      lambda_diff = max(abs(array_operator(lambda_i, lambda_i_upd, '-')));
      lambda_i = array_operator(0.05*lambda_i_upd, (1-0.05)*lambda_i, '+');
    }

    outerdiff = max(max(w_diff), max(lambda_diff))
    iter = iter + 1;

    # if(length(zeta_mult) > 0){
    #   if(length(zeta_mult) != length(zeta_intervals)){
    #     stop(
    #       "The array zeta_mult must have the same size as the array zeta_intervals."
    #     )
    #   } else{
    #     iiter = 1
    #     zeta_mult2 = c(zeta_mult, 1)
    #     while(outerdiff > zeta_intervals[iiter] & iiter <= (length(zeta_intervals))){
    #       iiter = iiter+1
    #     }
    #     zeta=zeta_init*zeta_mult[iiter]
    #     print(zeta_mult[iiter])
    #   }
    # }


    if(outerdiff<500){
      zeta = zeta_init*10;
    }

    if(outerdiff<300){
      zeta = zeta_init*50;
    }

    if(outerdiff<100){
      zeta = zeta_init*100;
    }

    if(outerdiff<20){
      zeta = zeta_init*500;
    }

    if(outerdiff<10){
      zeta = zeta_init*500;
    }

    if(outerdiff<1){
      zeta = zeta_init*500;
    }
    print(outerdiff)
  }

  # Save and export
  M_sum = sumDims(M,1);

  return(list(w=w, W_is=W_is, B=B, r=r, q=q, pi_js=pi_js, pi_s=pi_s, lambda_ijs_is=lambda_ijs_is,
              y_bar=y_bar, lambda_is_i=lambda_is_i, H_r=H_r, H_f=H_f, U=U, U_i=U_i,
              lambda_i=lambda_i, L_j=L_j, H_bar=H_bar, H_bar_rest=H_bar_rest, Q_js=Q_js,
              Q_s=Q_s, Q=Q, L_j_eff=L_j_eff, M_sum=M_sum))
}




