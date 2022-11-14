%% This file writes the mod files with the given parameters

%%
habit = 0.65;       % Parameter of habit formation
bet = 0.995;        % Discount factor
delta = 0.02;       % Depreciation of physical capital
reta = 0.001832219; % Parameter of leisure in utility
psi_l = 0.455;      % Parameter related to labour supply (labour disutility is: -psi_l * l^(1+elas_l)  / (1+elas_l) 
elas_l = 0.429;     % Elasticity of labour utility psi_l = 4.55 and elas_l=0.429  
pi_ss = 1.02^0.25;  % Steady state level of inflation
xi_w = 21;          % Parameter governing nominal rigidities and elasticities of substitution following CEE, 2005 xi_w = 21
phi_w = 0.64;       % Parameter governing nominal rigidities of wage setting

eta_h = 0.9;        % Measure of agents that are households
eta_e = 0.07;       % Measure of agents that are entrepreneurs
eta_b = 1-eta_h-eta_e; % Measure of agents that are bankers
xi_p = 6;           % Parameter governing nominal rigidities of the intermediate producers
phi_p = 0.6;        % Parameter of price updating probability of intermediate producers

theta_h = 0.639999; % Share of household labor in production
theta_k = 0.36;     % Share of capital in production
theta_e = (1.0-theta_h-theta_k)/2.0;       % Share of entrepreneurial labor in production
theta_b = 1.0-theta_h-theta_k-theta_e;     % Share of bank labor in production
bigR = 1.00;        % Overall return to the investment project 
delalpha = 0.35;    % Difference in probability of success between high and low effort
tau_b = 0.7;        % Survival probability of banker
tau_e = 0.7;        % Survival probability of entrepreneurs
lam_r = 0.9;        % AR(1) parameter in the Taylor rule, for short term interest rate smoothing
lam_pi = 1.19;       % Parameter on the inflation target in the Taylor rule 
lam_y = 0.176;      % Parameter of output stabilization in the Taylor rule
bby = 0.000;        % Additional capital endowment banks receive every period, as a percentage of output, set 0 by default.
nu = 1;             % Valuation mulitplier for exiting bankers consumption
rhoz = 0.95;        % AR(1) coefficient for the technology shock
rhomp = 0.0;        % AR(1) coefficient for the monetary policy
rhobk = 0.90;       % AR(1) coefficient of the bank capital shock
rhoqe = 0.9;        % AR(1) coefficient of the QE shock process
sigma_a = 0.5;	    % Parameter related to capital utilisation:

sigmaz = 0.0035;    % Standard deviation of the technology shock
sigmamp = 0.0016;   % Standard deviation of the monetary policy shock
sigmabk = 2.5;      % Standard deviation of the bank capital shock
alpha_ss = 0.993;   % Steady state probability of success given high effort that is also the bank success probability
mu_ss = 0.025;      % Steady state of monitoring effort by banks
omega = -5;         % Strength of response of monetary authority to deviations from credittogdp - initial guess to be updated by RP
varsigma = 0;       %0.01;   Strength of the endogenous link between credit to GDP devation from SS to banking sector riskiness - Parameter linking endogenous riskiness of the banking sector to total lending
epsb = 10;          % Linkage between shirking and premium paid  Parameter linking probability of success to monitoring

Btopbar =  11.503*0.9*0.35*1.05; %maximum private benefits from shirking it is not 0.1575 as in the table but Btopbar =delalpha*0.9*bigR; % Maximum level of private benefit

Blowbar = 0; 	    %minimum private benefit from shirking given monitoring -  Minimum level of private benefit
Chi =  15;            % 3.0792 Sensitivity of entrepreneurial private benefit to monitoring intensity - Parameter for endogenous private benefit
CBBSp = -0.81;       % Parameter governing the QE/Central bank balance sheet reaction
sigmaqe = 1;        % Standard deviation of QE shock