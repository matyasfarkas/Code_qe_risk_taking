// A Dynare version of Counter cyclical Monetary policy at the ZLB based on  
// Meh and Moran. 2010. "The Role of Bank Capital in the Propagation of Shocks" 
// JEDC 34: 555-576 
// Written by Matyas Farkas, 2019 
//--------------------------------------------------------------------- 
// 1. Variable declaration 
//--------------------------------------------------------------------- 
var  
ch    Ce    Cb    H     Y     I     bigN  bigA  smallb q     rk    infl  w_h   w_e   w_b   Rd    Ra    lam   mu     p     TL    mc    s     mgrowth expinfl smalld  gY   totC u    keff  K     Kb    Ke    ptilde nump   denp   wtilde numw   denw   lz     lmp    lbk    log_y  log_I  gammag alpha  rnot   CBBS lqe  util welf BVaR;
varexo  
z_shk  mp_shk bk_shk qe_shk macropru_shk; 
//--------------------------------------------------------------------- 
// 2. Parameter declaration and calibration 
//--------------------------------------------------------------------- 
parameters  
habit  bet      delta   reta     psi_l  elas_l chi1 chi2  pi_ss xi_w phi_w mark_w eta_h eta_e eta_b xi_p phi_p mark_p theta_k theta_h theta_e theta_b bigR delalpha tau_b tau_e lam_r  lam_pi lam_y bigtheta bby nu rhoz rhomp rhobk sigmaz sigmamp sigmabk  alpha_ss mu_ss omega varsigma epsb Btopbar Blowbar Chi CBBSp sigmaqe rhoqe ; 
sigmaqe = 0.00000000000000; 
rhoqe = 0.9000; 
habit = 0.6500;  bet = 0.9950;    delta = 0.0200;  reta = 0.00183221900000;
%%%%%%%%%%%%%%%%%%%%%%%%%
psi_l = 0.4550;
elas_l = 0.4290;
%%%%%%%%%%%%%%%%%%%%%%%%%
pi_ss = 1.0049629;
lam_r = 0.9000;
lam_pi = 1.8500;
lam_y = 0.0000;
Rd_temp = pi_ss/bet;
xi_w = 21.0;
phi_w = 0.6400;
xi_p = 6.0000;
phi_p = 0.6000;
mark_p  = xi_p/(xi_p-1);
mark_w = xi_w/(xi_w-1);
eta_h = 0.9000;            eta_e = 0.0700;           eta_b = 1-eta_h-eta_e;    theta_h = 0.63999900000000;   theta_k = 0.36000000000000;       theta_e = (1.0-theta_h-theta_k)/2.0;    theta_b = 1.0-theta_h-theta_k-theta_e;   delalpha = 0.35000000000000;
tau_b = 0.70000000000000; 
tau_e = 0.70000000000000;
bigR = 1.0200000000000; 
bby = 0.0000;
                                             nu = 1.00000000000000;
rhoz = 0.95000000000000;    sigmaz = 0.010000000000; sigma_a = 0.50000000000000;
rhomp = 0.00000000000000;  sigmamp = 0.010000000000; rhobk = 0.90000000000000;  sigmabk = 1.1500000000000; alpha_ss = 0.99300000000000;
Chi =  15.0000000000000;            epsb = 10.0000000000000;      Btopbar =  3.800; Blowbar = 0.0000; 
omega = -5.00000000000000;     varsigma = 0.01000000000000;  mu_ss = 0.02500000000000;
%%%%%%%%%%%%%%%%%%%%% 
alphatemp = alpha_ss; 
mutemp = mu_ss; 
smallbtemp = Btopbar*(1+Chi*mutemp)^(-epsb); 
kk1 = delta*theta_k/(alphatemp*bigR*(1/bet-1+delta));  
betatemp  = 1/bet; 
infl_temp = pi_ss; 
qtemp = kk1*(alphatemp*tau_b*mutemp*betatemp/delalpha-(1+mutemp+alphatemp*mutemp/(Rd_temp*delalpha))); 
qtemp = qtemp/(alphatemp*smallbtemp*kk1/(Rd_temp*delalpha)-alphatemp*bigR*kk1/Rd_temp-theta_e-theta_b-bby-alphatemp*tau_e*smallbtemp*kk1*betatemp/delalpha); 
Gtemp = 1+mutemp-(qtemp*alphatemp/Rd_temp)*(bigR-mutemp/(qtemp*delalpha)-smallbtemp/delalpha); 
IY = kk1/qtemp; 
KY = alphatemp*bigR*IY/delta;  
KeY = tau_e*alphatemp*smallbtemp*IY/delalpha; 
KbY = tau_b*alphatemp*mutemp*IY/(qtemp*delalpha); 
KhY = KY - KeY - KbY; 
CeY = (1-tau_e)*alphatemp*smallbtemp*IY*qtemp/delalpha; 
CbY = nu*(1-tau_b)*alphatemp*mutemp*IY/delalpha; 
ChY = 1 + bby - CeY - CbY  - (1+mutemp)*IY; 
essai = (1-bet*habit)*theta_h/((1-habit)*ChY*psi_l*mark_w); 
smallh = essai^(1/(1+elas_l)); 
Htemp = smallh*eta_h^(xi_w/(xi_w-1)); 
kk2 = (KY^(1/(1-theta_k))) * (eta_e^(theta_e/(1-theta_k)))* (eta_b^(theta_b/(1-theta_k)))* mark_p^(1/(theta_k-1)); 
Ktemp = Htemp^(theta_h/(1-theta_k))* kk2; 
Ytemp = (1/mark_p)* (Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b); 
bigtheta = (mark_p-1)*Ytemp; 
Itemp = IY*Ytemp; 
rktemp = (1/mark_p)*theta_k*(Ktemp^(theta_k-1))*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b); 
%%%%%%%%%%%% 
chi1 = rktemp; 
chi2 = sigma_a*chi1;   
%%%%%%%%%%%% 
w_htemp = (1/mark_p)*theta_h*(Ktemp^theta_k)*(Htemp^(theta_h-1))*(eta_e^theta_e)*(eta_b^theta_b); 
w_etemp = (1/mark_p)*theta_e*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^(theta_e-1))*(eta_b^theta_b); 
w_btemp = (1/mark_p)*theta_b*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^(theta_b-1)); 
Ketemp = KeY*Ytemp; 
Kbtemp = KbY*Ytemp; 
Khtemp = KhY*Ytemp; 
Cetemp = CeY*Ytemp; 
Cbtemp = CbY*Ytemp; 
chtemp = ChY*Ytemp/eta_h; 
bigAtemp = eta_b*w_btemp+(rktemp+qtemp*(1-delta))*Kbtemp+bby*Ytemp; 
bigNtemp = eta_e*w_etemp+(rktemp+qtemp*(1-delta))*Ketemp; 
totCtemp = Cetemp+Cbtemp+eta_h*chtemp; 
lamtemp = (1-bet*habit)/((1-habit)*(chtemp)); 
smalldtemp = alphatemp*qtemp*( bigR - smallbtemp/delalpha - mutemp/(qtemp*delalpha))*Itemp/(Rd_temp*eta_b); 
ptemp = (Rd_temp-1)*lamtemp/ ( (Rd_temp-1)*lamtemp*eta_b*smalldtemp + eta_h*reta); 
mctemp = reta*ptemp/((Rd_temp-1)*lamtemp); 
Ratemp = alphatemp*qtemp*(mutemp/(qtemp*delalpha))*Itemp/bigAtemp; 
TLtemp = Itemp-bigNtemp; 
wtildetemp = w_htemp*infl_temp*eta_h^(-1/(1-xi_w)); 
numwtemp =  ( ( (w_htemp*infl_temp)^(xi_w*(1+elas_l)) * Htemp^(1+elas_l) )/(1-bet*phi_w) )^(1/(1+elas_l*xi_w)); 
denwtemp = ( (w_htemp^(xi_w)*infl_temp^(xi_w-1)*Htemp*lamtemp)/(1-bet*phi_w) )^(1/(1+elas_l*xi_w)); 
wtildetemp = (mark_w*psi_l)^(1/(1+xi_w*elas_l))*numwtemp/denwtemp; 
stemp = 1/mark_p; 
mgrowthtemp = infl_temp; 
expinfltemp = infl_temp; 
gYtemp = 1; 
utemp = 1; 
gammagtemp = (Itemp-bigNtemp )/bigAtemp; 
%%%%%%%%%%%%%%%//--------------------------------------------------------------------- 
// 3. Model declaration 
//--------------------------------------------------------------------- 
model;   
1/(ch - habit*ch(-1)) - bet*habit*1/(ch(+1) - habit*ch) = lam; 
reta/(mc/p) = (Rd-1)*lam; 
rk = chi1 + chi2*(u - 1); 
q*lam = bet* (lam(+1)* ( q(+1)*(1-delta) + rk(+1)*u- ( chi1*(u-1)+0.5*chi2*(u(+1)-1)^2 ) )); 
lam = bet*lam(+1)*Rd(+1)/infl(+1); 
infl = p*mgrowth/p(-1); 
wtilde = (xi_w *psi_l/(xi_w-1))^(1/(1+elas_l*xi_w)) * numw /  denw; 
numw =  ( bet*phi_w*(numw(+1)^(1+elas_l*xi_w)) + (H^(1+elas_l)) *(w_h*infl)^(xi_w*(1+elas_l)) )^(1/(1+elas_l*xi_w)); 
denw =  ( bet*phi_w*(denw(+1)^(1+elas_l*xi_w)) + lam*H*(w_h^xi_w)*(infl^(xi_w-1)) )^(1/(1+elas_l*xi_w)) ; 
infl^(1-xi_w) *w_h^(1-xi_w) = eta_h*(1-phi_w)*wtilde^(1-xi_w) + phi_w* infl(-1)^(1-xi_w) *w_h(-1)^(1-xi_w); 
ptilde = xi_p * nump/ ((xi_p - 1) * denp ); 
nump = bet*phi_p*nump(+1) + lam*s*Y*infl^xi_p; 
denp = bet*phi_p*denp(+1) + lam*Y*infl^(xi_p-1); 
infl^(1-xi_p) = (1-phi_p)*ptilde^(1-xi_p) + phi_p*infl(-1)^(1-xi_p); 
rk = s * exp(lz)* theta_k * keff^(theta_k-1) * H^theta_h * eta_e^theta_e *eta_b^(theta_b); 
w_h = s * exp(lz)* theta_h * keff^(theta_k) * H^(theta_h-1) * eta_e^theta_e * eta_b^(theta_b); 
w_e = s * exp(lz)* theta_e * keff^(theta_k) * H^(theta_h) * eta_e^(theta_e-1) * eta_b^(theta_b); 
w_b = s * exp(lz)* theta_b * keff^(theta_k) * H^(theta_h) * eta_e^(theta_e) * eta_b^(theta_b-1); 
bigA*Ra = alpha * mu*I/delalpha; 
Rd*(1-eta_h*mc)/p = q*alpha*( bigR-smallb/delalpha - mu/(q*delalpha))*I; 
Kb = tau_b*alpha*mu*I/(q*delalpha); 
Ke = tau_e*alpha*smallb*I/delalpha; 
bigA = exp(lbk)*( rk + q*(1-delta))*Kb(-1) + eta_b*w_b; 
bigN = ( rk + q*(1-delta))*Ke(-1) + eta_e*w_e; 
Cb = (1-tau_b)*alpha*mu*I/delalpha; 
Ce = (1-tau_e)*alpha*q*smallb*I/delalpha; 
eta_h *ch + Ce + Cb + I + I*mu = Y; 
K = (1-delta)*K(-1) + alpha*bigR*I; 
log(rnot/steady_state(rnot)) = lam_r*log( rnot(-1)/steady_state(rnot) ) + (1-lam_r) * ( lam_pi*log(infl/steady_state(infl)) + lam_y*log(Y/steady_state(Y)) ) + lmp; 
Rd = rnot; 
Y = exp(lz)* (keff^theta_k)*(H^theta_h)*(eta_e^theta_e)*(eta_b^theta_b) - bigtheta; 
eta_b* smalld = (1-eta_h*mc)/p; 
keff = u*( K(-1) - Ke(-1) - Kb(-1) ) + Ke(-1) + Kb(-1); 
expinfl = infl(+1); 
gY = Y/Y(-1); 
TL = I - bigN; 
totC = eta_h *ch + Ce + Cb; 
lz = rhoz*lz(-1) -z_shk;  %+z_shk
lmp = rhomp*lmp(-1)+ mp_shk; 
lbk = rhobk*lbk(-1) -bk_shk; %+bk_shk
log_y = log(Y); 
log_I = log(I); 
smallb = Btopbar *(1+Chi*mu)^(-epsb); 
alpha = alpha_ss;
//alpha = 1 - (1-alpha_ss) - ((I-bigN)/Y-(steady_state(I)-steady_state(bigN))/steady_state(Y))*varsigma ;  
gammag = steady_state(gammag) + omega*((I-bigN-CBBS*Y)/Y - (steady_state(I)-steady_state(bigN))/steady_state(Y)) - macropru_shk;  
I =gammag * bigA + bigN + CBBS*Y; 
gammag = 1 + (q*Ra/Rd)*(delalpha*bigR-smallb-mu/q)/mu - Ra*delalpha/alpha;  
CBBS = lqe/steady_state(Y); 
lqe = rhoqe*lqe(-1)+ qe_shk; 
util = log( ch - habit*ch(-1) ) -psi_l * (1-bet*habit)*theta_h/((1-habit)*ch*psi_l*mark_w)  / (1+elas_l) + reta*log(mc/p); 
welf = util + bet*util(+1); 
BVaR = mu/(q*(exp(0.2345)-1));
end; 
//--------------------------------------------------------------------- 
// 4. Initial values and steady state 
//--------------------------------------------------------------------- 
resid; 
steady;   
check;                     
//--------------------------------------------------------------------- 
// 5. Shock declaration   
//                        
//--------------------------------------------------------------------- 
shocks; 
var z_shk = sigmaz^2; 
var mp_shk = sigmamp^2; 
var bk_shk = sigmabk^2; 
var qe_shk = sigmaqe^2; 
var macropru_shk = 0.01^2;
end; 
optim_weights; 
welf 1; 
end; 
osr_params omega; 
osr(opt_algo = 9, nograph); 
stoch_simul ch    Ce    Cb    H     Y     I     bigN  bigA  smallb q     rk    infl  w_h   w_e   w_b   Rd    Ra    lam   mu     p     TL    mc    s     mgrowth expinfl smalld  gY   totC u    keff  K     Kb    Ke    gammag alpha  rnot  util welf BVaR;