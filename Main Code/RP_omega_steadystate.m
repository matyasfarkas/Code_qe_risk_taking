% computes the steady state of BankCapital.mod. 
% Taken from an example file by S. Adjemian on Dynare Forum
% January 2011, Kevin Moran

function [ys,check] = BankCapital_steadystate(ys,exe)
global M_

%% DO NOT CHANGE THIS PART.
%%
%% Here we load the values of the deep parameters in a loop.
%%
NumberOfParameters = M_.param_nbr;                            % Number of deep parameters.
for i = 1:NumberOfParameters                                  % Loop...
  paramname = deblank(M_.param_names(i,:));                   %    Get the name of parameter i. 
  eval([ paramname ' = M_.params(' int2str(i) ');']);         %    Get the value of parameter i.
end                                                           % End of the loop.  
check = 0;
%%
%% END OF THE FIRST MODEL INDEPENDENT BLOCK.


%% THIS BLOCK IS MODEL SPECIFIC.
%%
%% Here the user has to define the steady state.
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% steady-state computations:

alpha = alpha_ss;

smallb_ss = Btopbar * (1+Chi*  mu_ss )^(-epsb);

mu_ss  =  mu_ss;

infl_ss = pi_ss;

Rd_ss = pi_ss/bet;

kk1 = delta*theta_k/(alpha_ss*bigR*(1/bet-1+delta)); 
betatemp  = 1/bet;

q_ss = kk1*(alpha_ss*tau_b*mu_ss*betatemp/delalpha-(1+mu_ss+alpha_ss*mu_ss/(Rd_ss*delalpha)));
q_ss = q_ss/(alpha_ss*smallb_ss*kk1/(Rd_ss*delalpha)-alpha_ss*bigR*kk1/Rd_ss-theta_e-theta_b-bby-alpha_ss*tau_e*smallb_ss*kk1*betatemp/delalpha);

%G_ss = 1+mu_ss-(q_ss*alpha_ss/Rd_ss)*(bigR-mu_ss/(q_ss*delalpha)-smallb_ss/delalpha);
   
IY = kk1/q_ss;

KY = alpha_ss*bigR*IY/delta; 

KeY = tau_e*alpha_ss*smallb_ss*IY/delalpha;
    
KbY = tau_b*alpha_ss*mu_ss*IY/(q_ss*delalpha);

KhY = KY - KeY - KbY;

CeY = (1-tau_e)*alpha_ss*smallb_ss*IY*q_ss/delalpha;

CbY = nu*(1-tau_b)*alpha_ss*mu_ss*IY/delalpha;

ChY = 1 + bby - CeY - CbY  - (1+mu_ss)*IY;

essai = (1-bet*habit)*theta_h/((1-habit)*ChY*psi_l*mark_w);

smallh = essai^(1/(1+elas_l));

H_ss = smallh*eta_h^(xi_w/(xi_w-1));

kk2 = (KY^(1/(1-theta_k))) * (eta_e^(theta_e/(1-theta_k)))* (eta_b^(theta_b/(1-theta_k)))* mark_p^(1/(theta_k-1));

K_ss = H_ss^(theta_h/(1-theta_k))* kk2;

Y_ss = (1/mark_p)* (K_ss^theta_k)*(H_ss^theta_h)*(eta_e^theta_e)*(eta_b^theta_b);

I_ss = IY*Y_ss;

rk_ss = (1/mark_p)*theta_k*(K_ss^(theta_k-1))*(H_ss^theta_h)*(eta_e^theta_e)*(eta_b^theta_b);

wh_ss = (1/mark_p)*theta_h*(K_ss^theta_k)*(H_ss^(theta_h-1))*(eta_e^theta_e)*(eta_b^theta_b);
                                 
we_ss = (1/mark_p)*theta_e*(K_ss^theta_k)*(H_ss^theta_h)*(eta_e^(theta_e-1))*(eta_b^theta_b);

wb_ss = (1/mark_p)*theta_b*(K_ss^theta_k)*(H_ss^theta_h)*(eta_e^theta_e)*(eta_b^(theta_b-1));

Ke_ss = KeY*Y_ss;

Kb_ss = KbY*Y_ss;

Kh_ss = KhY*Y_ss;

Ce_ss = CeY*Y_ss;

Cb_ss = CbY*Y_ss;

ch_ss = ChY*Y_ss/eta_h;

bigA_ss = eta_b*wb_ss+(rk_ss+q_ss*(1-delta))*Kb_ss+bby*Y_ss;
			
bigN_ss = eta_e*we_ss+(rk_ss+q_ss*(1-delta))*Ke_ss;

totC_ss = Ce_ss+Cb_ss+eta_h*ch_ss;

lam_ss = (1-bet*habit)/((1-habit)*(ch_ss));

smalld_ss = alpha_ss*q_ss*( bigR - smallb_ss/delalpha - mu_ss/(q_ss*delalpha))*I_ss/(Rd_ss*eta_b);
   
p_ss = (Rd_ss-1)*lam_ss/ ( (Rd_ss-1)*lam_ss*eta_b*smalld_ss + eta_h*reta);

mc_ss = reta*p_ss/((Rd_ss-1)*lam_ss);

Ra_ss = alpha_ss*q_ss*(mu_ss/(q_ss*delalpha))*I_ss/bigA_ss;

TL_ss = I_ss-bigN_ss;

%CA_ss = bigA_ss/((1+mu_ss)*I_ss-bigN_ss);
gammag_ss  = (I_ss-bigN_ss )/bigA_ss;

numw_ss =  ( ( (wh_ss*infl_ss)^(xi_w*(1+elas_l)) * H_ss^(1+elas_l) )/(1-bet*phi_w) )^(1/(1+elas_l*xi_w));

denw_ss = ( (wh_ss^(xi_w)*infl_ss^(xi_w-1)*H_ss*lam_ss)/(1-bet*phi_w) )^(1/(1+elas_l*xi_w));

wtilde_ss = (mark_w*psi_l)^(1/(1+xi_w*elas_l))*numw_ss/denw_ss;

s_ss = 1/mark_p;

mgrowth_ss = infl_ss;

expinfl_ss = infl_ss;

gY_ss = 1;

u_ss = 1;

gammag_ss  = (I_ss-bigN_ss )/bigA_ss;
% Now ready to assign these "*_ss values to all variables
s = s_ss; 
u = u_ss;
infl = infl_ss;
expinfl = expinfl_ss;
mgrowth = mgrowth_ss;
gY = gY_ss;
Rd = Rd_ss; 
q = q_ss;
% = G_ss;
H = H_ss;
K = K_ss;
Y = Y_ss;
I = I_ss;          
rk = rk_ss;
w_h = wh_ss;
w_e = we_ss;
w_b = wb_ss;
Ke = Ke_ss;
Kb = Kb_ss;
Ce = Ce_ss;
Cb = Cb_ss;
ch = ch_ss;
totC = totC_ss;
lam = lam_ss;
bigA = bigA_ss;			
bigN = bigN_ss;
smalld = smalld_ss;   
p = p_ss;
mc = mc_ss;
Ra = Ra_ss;
TL = TL_ss; 
numw = numw_ss;
denw = denw_ss;
wtilde = wtilde_ss;
nump = (lam_ss * Y_ss * s_ss * infl_ss^xi_p)/(1-bet*phi_p);
denp = (lam_ss * Y_ss * infl_ss^(xi_p-1))/(1-bet*phi_p);
ptilde = xi_p * nump/ ((xi_p - 1) * denp );
lz = 0.0;
lmp = 0.0;
lbk = 0.0;
keff = K_ss;
%CA = CA_ss;
log_y = log(Y);
log_I = log(I);
gammag=gammag_ss;
smallb=smallb_ss;
mu=mu_ss;
alpha = alpha_ss;  % endogenous probability of default
rnot = Rd_ss;
CBBS = 0;
lqe = 0;
bdp =mu_ss/q_ss;
%inflexp10 = infl_ss;
inflexp1 = infl_ss;
inflexp2 = inflexp1;
inflexp3 = inflexp2;
inflexp4 = inflexp3;
inflexp5 = inflexp4;
inflexp6 = inflexp5;
inflexp7 = inflexp6;
inflexp8 = inflexp7;
inflexp9 = inflexp8;
inflexp10 = inflexp9;
inflexp11 = inflexp10;
inflexp12 = inflexp11;
inflexp13 = inflexp12;
inflexp14 = inflexp13;
inflexp15 = inflexp14;
inflexp16 = inflexp15;
inflexp17 = inflexp16;
inflexp18 = inflexp17;
inflexp19 = inflexp18;
inflexp20 = inflexp19;
inflexp21 = inflexp20;
inflexp22 = inflexp21;
inflexp23 = inflexp22;
inflexp24 = inflexp23;
inflexp25 = inflexp24;
inflexp26 = inflexp25;
inflexp27 = inflexp26;
inflexp28 = inflexp27;
inflexp29 = inflexp28;
inflexp30 = inflexp29;
inflexp31 = inflexp30;
inflexp32 = inflexp31;
inflexp33 = inflexp32;
inflexp34 = inflexp33;
inflexp35 = inflexp34;
inflexp36 = inflexp35;
inflexp37 = inflexp36;
inflexp38 = inflexp37;
inflexp39 = inflexp38;
inflexp40 = inflexp39;



Rexp1 = Rd_ss;
Rexp2 = Rexp1;
Rexp3 = Rexp2;
Rexp4 = Rexp3;
Rexp5 = Rexp4;
Rexp6 = Rexp5;
Rexp7 = Rexp6;
Rexp8 = Rexp7;
Rexp9 = Rexp8;
Rexp10 = Rexp9;
Rexp11 = Rexp10;
Rexp12 = Rexp11;
Rexp13 = Rexp12;
Rexp14 = Rexp13;
Rexp15 = Rexp14;
Rexp16 = Rexp15;
Rexp17 = Rexp16;
Rexp18 = Rexp17;
Rexp19 = Rexp18;
Rexp20 = Rexp19;
Rexp21 = Rexp20;
Rexp22 = Rexp21;
Rexp23 = Rexp22;
Rexp24 = Rexp23;
Rexp25 = Rexp24;
Rexp26 = Rexp25;
Rexp27 = Rexp26;
Rexp28 = Rexp27;
Rexp29 = Rexp28;
Rexp30 = Rexp29;
Rexp31 = Rexp30;
Rexp32 = Rexp31;
Rexp33 = Rexp32;
Rexp34 = Rexp33;
Rexp35 = Rexp34;
Rexp36 = Rexp35;
Rexp37 = Rexp36;
Rexp38 = Rexp37;
Rexp39 = Rexp38;
Rexp40 = Rexp39;


util_ss = log( ch_ss - habit*ch_ss ) -psi_l * (1-bet*habit)*theta_h/((1-habit)*ch_ss*psi_l*mark_w)  / (1+elas_l) + reta*log(mc_ss/p_ss);
                
welf_ss = util_ss + bet*util_ss;

util = util_ss;
welf = welf_ss;


%%
%% END OF THE MODEL SPECIFIC BLOCK.


%% DO NOT CHANGE THIS PART.
%%
%% Here we define the steady state values of the endogenous variables of
%% the model.
%%
NumberOfEndogenousVariables = M_.endo_nbr;                    % Number of endogenous variables.
ys = zeros(NumberOfEndogenousVariables,1);                    % Initialization of ys (steady state).
for i = 1:NumberOfEndogenousVariables                         % Loop...
  varname = deblank(M_.endo_names(i,:));                      %    Get the name of endogenous variable i.                     
  eval(['ys(' int2str(i) ') = ' varname ';']);                %    Get the steady state value of this variable.
end                                                           % End of the loop.
%%
%% END OF THE SECOND MODEL INDEPENDENT BLOCK.