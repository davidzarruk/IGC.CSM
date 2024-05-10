# Load required libraries
library(shiny)
library(DT)  # For rendering the table

source("R/array_operator.R")
source("R/av_income_simple.R")
source("R/commuting_matrix.R")
source("R/density_development.R")
source("R/living_amenities_simple.R")
source("R/productivity.R")
source("R/solveModel.R")
source("R/sumDims.R")
source("R/sumDims2.R")


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
                          maxiter=1000){
  
  # Formatting of input data
  if(is.data.frame(L_i)){
    L_i = array(unlist(L_i), dim(L_i))
  } else if(is.null(dim(L_i))){
    L_i = array(L_i, dim=c(N,1))
  }
  
  if(is.data.frame(L_j)){
    L_j = array(unlist(L_j), dim(L_j))
  } else if(is.null(dim(L_j))){
    L_j = array(L_j, dim=c(N,1))
  }
  if(is.data.frame(K)){
    K = array(unlist(K), dim(K))  
  } else if(is.null(dim(K))){
    K = array(K, dim=c(N,1))
  }
  if(is.data.frame(Q)){
    Q = array(unlist(Q), dim(Q))
  } else if(is.null(dim(Q))){
    Q = array(Q, dim=c(N,1))
  }
  t_ij = array(unlist(t_ij), dim(t_ij))  
  
  # Normalize L_i to have the same size as L_j
  L_i=L_i*sum(L_j)/sum(L_i)
  
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

wages_inversion = function(N,
                           w_init,
                           theta,
                           tau,
                           L_i,
                           L_j,
                           nu_init=0.05,
                           tol=10^-10,
                           maxiter=10000){
  
  # Settings
  outerdiff = Inf
  w = w_init
  iter = 0
  nu = nu_init
  
  message("Inverting wages...\n")
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
    
    if(iter %% 10 == 0){
      message(paste0("Iteration: ", iter, ", error: ", round(outerdiff, 10), ".\n"))
    }
  }
  if(outerdiff<=tol){
    message(paste0("Converged after ", iter, " iterations. Error=", round(outerdiff, 10), ".\n"))
  } else{
    message(paste0("Reached maximum number of iterations (", iter, "). Error=", round(outerdiff, 10), ".\n"))
  }
  
  return(list(w=w, w_tr=w_tr, W_i=W_i, lambda_ij_i=lambda_ij_i))
}


# Define the UI for the Shiny app
ui <- fluidPage(
  shinyjs::useShinyjs(),
  titlePanel("Invert wages:"),
  sidebarLayout(
    sidebarPanel(
      # File input for the CSV dataset
      fileInput("data_locations", "Choose Data Locations (CSV):", accept = ".csv"),
      fileInput("data_times", "Choose Matrix Travel Times (CSV)", accept = ".csv"),
      numericInput("input_maxiter", "Number of iterations:", value = 100, min = 1, max = 10000, step = 1),
      
      # Execute button
      actionButton("execute", "Execute")
    ),
    mainPanel(
      # Output to display the result
      textOutput("text"),
      verbatimTextOutput("result"),
      
      # Output to display the selected dataset (optional)
      DTOutput("table")
    )
  )
)


# Define the server logic
server <- function(input, output, session) {
  # Increase the maximum upload size limit to 50 MB (52428800 bytes)
  options(shiny.maxRequestSize = 50*1024^2)
  
  # Read the selected CSV file
  dataset_data_locations <- reactive({
    req(input$data_locations)
    read.csv(input$data_locations$datapath, stringsAsFactors = FALSE)
  })

  # Read the selected CSV file
  dataset_data_times <- reactive({
    req(input$data_times)
    read.csv(input$data_times$datapath, stringsAsFactors = FALSE)
  })
  
  # React to the "Execute" button
  observeEvent(input$execute, {
    # Get the input values
    data_locations <- dataset_data_locations()
    data_times <- dataset_data_times()
    
    L_j = data_locations$t_w_vodacom # cantidad de trabajadores en cada location j
    L_i = data_locations$t_r_vodacom # cantidad de habitantes en cada location j
    L_i = L_i*sum(L_j)/sum(L_i)
    K = data_locations$SAL_area # tamaño del lugar
    Q = data_locations$price_m2
    N = length(L_i)
    t_ij = as.matrix(data_times[,2:(N+1)], dim=c(N,N))
    t_ij[864,] = t_ij[863,]
    t_ij[,864] = t_ij[,863]
    t_ij[,714] = t_ij[,713]
    t_ij[714,] = t_ij[713,]
    t_ij[165,] = t_ij[164,]
    t_ij[,165] = t_ij[,164]
    t_ij[,138] = t_ij[,137]
    t_ij[138,] = t_ij[137,]
    
    withCallingHandlers({
      shinyjs::html("text", "")
      # Call your function and capture its output
      inversionModel(N=N,
                     L_i=L_i,
                     L_j=L_j,
                     Q=Q,
                     K=K,
                     t_ij=t_ij,
                     maxiter=input$input_maxiter)
    },
    message = function(m) {
      shinyjs::html(id = "text", html = paste0(m$message, '<br>'), add = TRUE)
    })

  })
}

# Run the Shiny app
shinyApp(ui, server)