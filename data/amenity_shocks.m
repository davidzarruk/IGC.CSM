%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   m1_baseline_model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
name = 'metro' %% {metro, monorail, etrain, all, mon_etrain}
option = 2 %% 1 same agglomeration forces 2 different agglomeration forces

%% Preliminaries
if strcmp(getenv('USER'),'wb518521')==0
    cd("C:\Users\wb518521\Dropbox\Cairo\model\simulations_sector")
else
    cd('');
end

addpath(genpath('depends'));
rng(2020);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (0) Number of groups and sectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of shriakas
N = 875;
S = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import distance matrix 
dist_mat = readmatrix("C:\Users\wb518521\Dropbox\Cairo\temp\distance_matrix.txt");
dist_mat = dist_mat(:,2:876);

% Calculate travel time. Assuming speed of 15 km/h
speed = 10;
time_mat = (dist_mat ./ speed)*60;

% Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Workers and population
L_bar  = 1;
Data = importdata('C:\Users\wb518521\Dropbox\Cairo\temp\data_shriaka.csv');
Jobs = importdata('C:\Users\wb518521\Dropbox\Cairo\temp\employment_shriaka.csv');
Resi = importdata('C:\Users\wb518521\Dropbox\Cairo\temp\residents_shriaka.csv');
pop_t = Data.data(:,strcmp(Data.colheaders,'pop'));
workers_t = Jobs.data(:,strcmp(Jobs.colheaders,'total_empl_shriaka'));
workers_s = Jobs.data(:,3:5)*L_bar; 
pop = pop_t/sum(pop_t,1);
workers = workers_t./sum(workers_t,1);
lambda_is_i = Resi.data(:,2:4);
res_agg_s = sum(pop.*lambda_is_i,1);
wor_agg_s = sum(workers_s,1);
workers_s = (workers_s.*res_agg_s)./wor_agg_s;
%workers_s = workers;
%res_agg_s = 1;
%lambda_is_i = ones(N,1);
FAR = Data.data(:,strcmp(Data.colheaders,'far'));
area = Data.data(:,strcmp(Data.colheaders,'Area'))./1000000;
%area = Data(:,7)*0.00001;

%% New and Old cities
NC = importdata('C:\Users\wb518521\Dropbox\Cairo\temp\new_cities.csv');
new = NC.data(:,2);
OC = importdata('C:\Users\wb518521\Dropbox\Cairo\temp\old_cities.csv');
old = OC.data(:,2);
new_old = new+old;
AmenityShock = [old new new new_old new];

% distance to transit intervenstions
Distance = importdata("C:\Users\wb518521\Dropbox\Cairo\temp\distance_transit.txt");
distance_metro = Distance(:,2);
distance_hsr = Distance(:,3);
distance_monorail = Distance(:,4);
distance_etrain = Distance(:,5);

% Iceberg commuting cost
lambda = 0.01;
tau = exp(lambda * time_mat);
discount1 = 1-((exp(0.2*(distance_metro/150))).^(-1))*0.1;
discount_tr = discount1';
tau_shock_metro = max((tau.*discount1).*discount_tr,1);
tau_hat_metro = tau_shock_metro./tau;

discount2 = 1-((exp(0.25*(distance_monorail/150))).^(-1))*0.12;
discount_tr = discount2';
tau_shock_monorail = max((tau.*discount2).*discount_tr,1);
tau_hat_monorail = tau_shock_monorail./tau;

discount3 = 1-((exp(0.25*(distance_etrain/150))).^(-1))*0.195;
discount_tr = discount3';
tau_shock_etrain = max((tau.*discount3).*discount_tr,1);
tau_hat_etrain = tau_shock_etrain./tau;

discount4 = discount1.*discount2.*discount3;
discount_tr = discount4';
tau_shock_all = max((tau.*discount4).*discount_tr,1);
tau_hat_all = tau_shock_all./tau;

discount5 = discount2.*discount3;
discount_tr = discount5';
tau_shock_mon_etrain = max((tau.*discount5).*discount_tr,1);
tau_hat_mon_etrain = tau_shock_mon_etrain./tau;

%% Parameters
alpha1 = 0.7;
u = ones(N,1);
u_init = u;
% Lambda
A = ones(N,S);
B = ones(N,S);
theta1 = 4;
eta1 = 1.5;
kappa1 = 2;
xi1 = 1.8;
% Floor space supply
delta1 = 0.65*ones(N,1);
A_tilde = FAR;
T = ones(N,1);
% Number of firms
beta1 = 0.7*ones(1,1,S);
beta_tilde = (beta1.^(-beta1)).*((1-beta1).^(-(1-beta1)));
F = 1;
if option==1
    sigma1 = 6*ones(1,1,S);
end
if option==2
    sigma1 = 6*ones(1,1,S);
    sigma1(:,:,1) = 4;
    sigma1(:,:,2) = 4;
end
theta_set = [4 7 10];
percent_set = linspace(0.01,0.20,20);
% Residents in each location
lambda_i = pop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Solve Models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Note: procedure is: (i) solve model in GE to get w, r, and Lr 
%%
% Setup
%%
% Basic Settings
tol = 1e-10;
zeta = 0.0001;
A_init = A;
B_init = B;
w_init = ones(N,S);
r_init = ones(N,1);
q_init = ones(N,1);
vartheta_init = 0.5*ones(N,1);
H_bar = FAR.*area;
H_bar_rest = 12;

%% No migration
U_bl =[];
U = [];
Uhat_m = [];
lambda_hat_m = [];
L_hat_m = [];
Uhat_m_agg =[];
Qhat_m_agg =[];
Uhat_nm = [];
tau_init = tau;
H_bar_init = H_bar;

%% Test trade costs
tau = tau_init;
theta1=theta_set(2);
endo_Lr = 1;
param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);
data     = struct('tau',tau,'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);
settings = struct('tol',tol,'zeta',zeta,'w_init',w_init,'r_init',r_init,'q_init',q_init,'vartheta_init',vartheta_init,'endo_Lr',endo_Lr);
inversion_m_bl  = inversionModel_Eff(param,data,settings);
data.A = inversion_m_bl.A;
data.u = inversion_m_bl.u;
data.u_init = inversion_m_bl.u;
data.B = inversion_m_bl.B;
settings.w_init = inversion_m_bl.w;
results_m_bl  = solveModel1_Eff(param,data,settings);
namelist = ["metro"  "monorail" "etrain" "all" "mon_etrain"];
tau_shock = cat(3,tau_shock_metro,tau_shock_monorail,tau_shock_etrain,tau_shock_all,tau_shock_mon_etrain);
a = 1-cellfun('isempty',strfind(namelist,name));
index = find(a);

% Transit Intervention
A = inversion_m_bl.A;
u = inversion_m_bl.u;
B = inversion_m_bl.B;
w_init = inversion_m_bl.w;
endo_Lr = 1;
param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);
data     = struct('tau',tau_shock(:,:,index),'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);
settings = struct('tol',tol,'zeta',zeta,'w_init',w_init,'r_init',r_init,'q_init',q_init,'vartheta_init',vartheta_init,'endo_Lr',endo_Lr);
results_m_shock  = solveModel1_Eff(param,data,settings);
U_shock = results_m_shock.U_i./results_m_bl.U_i;
lambda_shock = results_m_shock.lambda_i./results_m_bl.lambda_i;
L_shock = results_m_shock.L_j./results_m_bl.L_j;
U_shock_agg = results_m_shock.U/results_m_bl.U;
Q_shock_agg = results_m_shock.Q_s./results_m_bl.Q_s;
Uhat_m = [Uhat_m U_shock];
lambda_hat_m = [lambda_hat_m lambda_shock];
L_hat_m = [L_hat_m L_shock];
Uhat_m_agg = [Uhat_m_agg U_shock_agg];
Qhat_m_agg =[Qhat_m_agg Q_shock_agg];

%% Amenity shock
A = inversion_m_bl.A;
u = inversion_m_bl.u.*(1+0.9*AmenityShock(:,index));
B = inversion_m_bl.B;
w_init = inversion_m_bl.w;
endo_Lr = 1;
param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);
data     = struct('tau',tau,'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);
settings = struct('tol',tol,'zeta',zeta,'w_init',w_init,'r_init',r_init,'q_init',q_init,'vartheta_init',vartheta_init,'endo_Lr',endo_Lr);
results_m_shock  = solveModel1_Eff(param,data,settings);
U_shock = results_m_shock.U_i./results_m_bl.U_i;
lambda_shock = results_m_shock.lambda_i./results_m_bl.lambda_i;
L_shock = results_m_shock.L_j./results_m_bl.L_j;
U_shock_agg = results_m_shock.U/results_m_bl.U;
Q_shock_agg = results_m_shock.Q_s./results_m_bl.Q_s;
Uhat_m = [Uhat_m U_shock];
lambda_hat_m = [lambda_hat_m lambda_shock];
L_hat_m = [L_hat_m L_shock];
Uhat_m_agg = [Uhat_m_agg U_shock_agg];
Qhat_m_agg =[Qhat_m_agg Q_shock_agg];

%% Amenity shock + Transit intervention
A = inversion_m_bl.A;
u = inversion_m_bl.u.*(1+0.9*AmenityShock(:,index));
B = inversion_m_bl.B;
w_init = inversion_m_bl.w;
endo_Lr = 1;
param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);
data     = struct('tau',tau_shock(:,:,index),'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);
settings = struct('tol',tol,'zeta',zeta,'w_init',w_init,'r_init',r_init,'q_init',q_init,'vartheta_init',vartheta_init,'endo_Lr',endo_Lr);
results_m_shock  = solveModel1_Eff(param,data,settings);
U_shock = results_m_shock.U_i./results_m_bl.U_i;
lambda_shock = results_m_shock.lambda_i./results_m_bl.lambda_i;
L_shock = results_m_shock.L_j./results_m_bl.L_j;
U_shock_agg = results_m_shock.U/results_m_bl.U;
Q_shock_agg = results_m_shock.Q_s./results_m_bl.Q_s;
Uhat_m = [Uhat_m U_shock];
lambda_hat_m = [lambda_hat_m lambda_shock];
L_hat_m = [L_hat_m L_shock];
Uhat_m_agg = [Uhat_m_agg U_shock_agg];
Qhat_m_agg =[Qhat_m_agg Q_shock_agg];

%% Changing the restriction by x% and changing amenities
per_rest_set = [0.5];
for i=1:1
    per_rest = per_rest_set(1,i);
    H_bar_rest = 12*(1+per_rest);
    u = inversion_m_bl.u;
    
    param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);
    data     = struct('tau',tau,'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);
    settings = struct('tol',tol,'zeta',zeta,'w_init',w_init,'r_init',r_init,'q_init',q_init,'vartheta_init',vartheta_init,'endo_Lr',endo_Lr);
    results_m_bl_shock  = solveModel1_Eff(param,data,settings);    
    U_transit_agg_shock = results_m_bl_shock.U/results_m_bl.U;
    Q_transit_agg_shock = results_m_bl_shock.Q_s./results_m_bl.Q_s;
 
    % Transit
    data     = struct('tau',tau_shock(:,:,index),'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);    
    results_m_transit  = solveModel1_Eff(param,data,settings);
    U_transit = results_m_transit.U_i./results_m_bl.U_i;
    lambda_transit = results_m_transit.lambda_i./results_m_bl.lambda_i;
    L_transit = results_m_transit.L_j./results_m_bl.L_j;
    U_transit_agg = (results_m_transit.U/results_m_bl.U)-U_transit_agg_shock;
    Q_transit_agg = (results_m_transit.Q_s./results_m_bl.Q_s)-Q_transit_agg_shock;
    Uhat_m = [Uhat_m U_transit];
    lambda_hat_m = [lambda_hat_m lambda_transit];
    L_hat_m = [L_hat_m L_transit];
    Uhat_m_agg = [Uhat_m_agg U_transit_agg];
    Qhat_m_agg = [Qhat_m_agg Q_transit_agg];
    
    % Amenity
    param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);    
    u = inversion_m_bl.u.*(1+0.9*AmenityShock(:,index)); 
    data     = struct('tau',tau,'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);    
    results_m_transit  = solveModel1_Eff(param,data,settings);
    U_transit = results_m_transit.U_i./results_m_bl.U_i;
    lambda_transit = results_m_transit.lambda_i./results_m_bl.lambda_i;
    L_transit = results_m_transit.L_j./results_m_bl.L_j;
    U_transit_agg = (results_m_transit.U/results_m_bl.U)-U_transit_agg_shock;
    Q_transit_agg = (results_m_transit.Q_s./results_m_bl.Q_s)-Q_transit_agg_shock;
    Uhat_m = [Uhat_m U_transit];
    lambda_hat_m = [lambda_hat_m lambda_transit];
    L_hat_m = [L_hat_m L_transit];
    Uhat_m_agg = [Uhat_m_agg U_transit_agg];
    Qhat_m_agg = [Qhat_m_agg Q_transit_agg];  
    
    % Amenity + Transit Intervention
    param    = v2struct(alpha1,beta1,beta_tilde,delta1,theta1,eta1,kappa1,xi1,sigma1,L_bar,H_bar,H_bar_rest,N,S);    
    u = inversion_m_bl.u.*(1+0.9*AmenityShock(:,index)); 
    data     = struct('tau',tau_shock(:,:,index),'A',A,'A_init',A_init,'u_init',u,'B',B,'B_init',B_init,'F',F,'A_tilde',A_tilde,'T',T,'lambda_i',lambda_i,'lambda_is_i',lambda_is_i,'L_j_data',workers_s);    
    results_m_transit  = solveModel1_Eff(param,data,settings);
    U_transit = results_m_transit.U_i./results_m_bl.U_i;
    lambda_transit = results_m_transit.lambda_i./results_m_bl.lambda_i;
    L_transit = results_m_transit.L_j./results_m_bl.L_j;
    U_transit_agg = (results_m_transit.U/results_m_bl.U)-U_transit_agg_shock;
    Q_transit_agg = (results_m_transit.Q_s./results_m_bl.Q_s)-Q_transit_agg_shock;
    Uhat_m = [Uhat_m U_transit];
    lambda_hat_m = [lambda_hat_m lambda_transit];
    L_hat_m = [L_hat_m L_transit];
    Uhat_m_agg = [Uhat_m_agg U_transit_agg];
    Qhat_m_agg = [Qhat_m_agg Q_transit_agg];      
end    
lambda_bl = results_m_bl.lambda_i.*results_m_bl.lambda_is_i;
lambda_shock = results_m_transit.lambda_i.*results_m_transit.lambda_is_i;
lambda_set = [lambda_bl lambda_shock];

option_str = string(option);
fileU = strcat('results/Uhat_transit_',name,option_str,'.csv');
fileL = strcat('results/L_hat_transit_',name,option_str,'.csv');
fileR = strcat('results/lambda_hat_transit_',name,option_str,'.csv');
fileUagg = strcat('results/Uhat_transit_agg_',name,option_str,'.csv');
fileQagg = strcat('results/Qhat_transit_agg_',name,option_str,'.csv');
fileLambda = strcat('results/Lambda_transit_',name,option_str,'.csv');

csvwrite(fileU,Uhat_m);
csvwrite(fileL,L_hat_m);
csvwrite(fileR,lambda_hat_m);
csvwrite(fileUagg,Uhat_m_agg);
csvwrite(fileQagg,Qhat_m_agg);
csvwrite(fileLambda,lambda_set);

