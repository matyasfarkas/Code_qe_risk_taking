%% Opening the non binding model file
 fid = fopen(['RP_omega.mod'],'w+');


   fprintf(fid,'// A Dynare version of Counter cyclical Monetary policy at the ZLB based on  \n');
   fprintf(fid,'// Meh and Moran. 2010. "The Role of Bank Capital in the Propagation of Shocks" \n');
   fprintf(fid,'// JEDC 34: 555-576 \n');
   fprintf(fid,'// Written by Matyas Farkas, 2019 \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'// 1. Variable declaration \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'var  \n');
   fprintf(fid,'ch    Ce    Cb    H     Y     I     bigN  bigA  smallb q     rk    infl  w_h   w_e   w_b   Rd    Ra    lam   mu     p     TL    mc    s     mgrowth expinfl smalld  gY   totC u    keff  K     Kb    Ke    ptilde nump   denp   wtilde numw   denw   lz     lmp    lbk    log_y  log_I  gammag alpha  rnot   CBBS lqe  bdp  inflexp1  inflexp2  inflexp3  inflexp4  inflexp5  inflexp6  inflexp7  inflexp8  inflexp9  inflexp10  inflexp11  inflexp12  inflexp13  inflexp14  inflexp15  inflexp16  inflexp17  inflexp18  inflexp19  inflexp20  inflexp21  inflexp22  inflexp23  inflexp24  inflexp25  inflexp26  inflexp27  inflexp28  inflexp29  inflexp30  inflexp31  inflexp32  inflexp33  inflexp34  inflexp35  inflexp36  inflexp37  inflexp38  inflexp39  inflexp40  Rexp1  Rexp2  Rexp3  Rexp4  Rexp5  Rexp6  Rexp7  Rexp8  Rexp9  Rexp10  Rexp11  Rexp12  Rexp13  Rexp14  Rexp15  Rexp16  Rexp17  Rexp18  Rexp19  Rexp20  Rexp21  Rexp22  Rexp23  Rexp24  Rexp25  Rexp26  Rexp27  Rexp28  Rexp29  Rexp30  Rexp31  Rexp32  Rexp33  Rexp34  Rexp35  Rexp36  Rexp37  Rexp38  Rexp39  Rexp40  util welf; \n');
   fprintf(fid,'varexo  \n');
   fprintf(fid,'z_shk  mp_shk bk_shk qe_shk ; \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'// 2. Parameter declaration and calibration \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'parameters  \n');
   fprintf(fid,'habit  bet      delta   reta     psi_l  elas_l chi1 chi2  pi_ss xi_w phi_w mark_w eta_h eta_e eta_b xi_p phi_p mark_p theta_k theta_h theta_e theta_b bigR delalpha tau_b tau_e lam_r  lam_pi lam_y bigtheta bby nu rhoz rhomp rhobk sigmaz sigmamp sigmabk  alpha_ss mu_ss omega varsigma epsb Btopbar Blowbar Chi CBBSp sigmaqe rhoqe ; \n');
  
%% Calibration
                fprintf(fid,'sigmaqe = %15.14f; \n',0);
                fprintf(fid,'rhoqe = %5.4f; \n',rhoqe);
                fprintf(fid,'habit = %5.4f;  %parameter governing habit formation (CEE 2005 use 0.65)\n', habit);
                fprintf(fid,'bet = %5.4f;    %discount factor\n', bet);
                fprintf(fid,'delta = %5.4f;  %depreciation of physical capital\n', delta);
                fprintf(fid,'reta = %15.14f;\n',reta);
                fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
                fprintf(fid,'% parameters related to labour supply\n');
                fprintf(fid,'% (labour disutility is: -psi_l * l^(1+elas_l)  / (1+elas_l) )\n');
                fprintf(fid,'% elas_l = 1 follows CEE(2005)\n');
                fprintf(fid,'% elas_l = 0.429 reproduces (for dynamics) our JEDC specification with log disutility\n');
                fprintf(fid,'% in either case, psi_l is set in order for work effort at ss to be 0.3\n');
                fprintf(fid,'% (this requires psi_l = 4.55 if elas_l = 0.429)\n');
                fprintf(fid,'%psi_l = 9.05;\n');
                fprintf(fid,'%elas_l = 1.0;\n');
                fprintf(fid,'psi_l = %5.4f;\n', psi_l);
                fprintf(fid,'elas_l = %5.4f;\n', elas_l);
                fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
                fprintf(fid,'%parameters of monetary policy\n');
                fprintf(fid,'pi_ss = %8.7f;\n', pi_ss);
                fprintf(fid,'lam_r = %5.4f;\n', lam_r);
                fprintf(fid,'lam_pi = %5.4f;\n', lam_pi);
                fprintf(fid,'lam_y = %5.4f;\n', lam_y);
                fprintf(fid,'Rd_temp = pi_ss/bet;\n');
                fprintf(fid,'% nominal rigidities and elasticities of substitution\n');
                fprintf(fid,'% (follows CEE, 2005)\n');
                fprintf(fid,'xi_w = %3.1f;\n', xi_w);
                fprintf(fid,'phi_w = %5.4f;\n', phi_w);
                fprintf(fid,'xi_p = %5.4f;\n',xi_p);
                fprintf(fid,'phi_p = %5.4f;\n',phi_p);
                fprintf(fid,'mark_p  = xi_p/(xi_p-1);\n');
                fprintf(fid,'mark_w = xi_w/(xi_w-1);\n');
                fprintf(fid,'eta_h = %5.4f;            % measure of agents that are households\n', eta_h);
                fprintf(fid,'eta_e = %5.4f;           % measure of agents that are entrepreneurs\n', eta_e);
                fprintf(fid,'eta_b = 1-eta_h-eta_e;    % measure of agents that are bankers\n');
                fprintf(fid,'theta_h = %15.14f;   % share of household labor in production\n', theta_h);
                fprintf(fid,'theta_k = %15.14f;       % share of capital in production\n', theta_k);
                fprintf(fid,'theta_e = (1.0-theta_h-theta_k)/2.0;    % share of entrepreneurial labor in production\n');
                fprintf(fid,'theta_b = 1.0-theta_h-theta_k-theta_e;   % share of bank labor in production\n');
                fprintf(fid,'% sector producing physical capital\n');
                fprintf(fid,'delalpha = %15.14f;\n',delalpha);
                fprintf(fid,'tau_b = %15.14f; \n',tau_b);
                fprintf(fid,'tau_e = %15.14f;\n', tau_e);
                fprintf(fid,'bigR = %15.14f; \n', bigR);
                fprintf(fid,'% In some experiments, we want to leave mu_ss unchanged at 0.025 (to compute s.s.)\n');
                fprintf(fid,'% but have mu = 0.0 for the linearization (to produce "WITH BANK CAPITAl CHANNEL" responses\n');
                fprintf(fid,'% and "NO BANK CAPITAl CHANNEL" responses)\n');
                fprintf(fid,'% in order to do that, need to disable lines 81-84 in the "steady_.m" program\n');
                fprintf(fid,'% of your dynare build\n');
                fprintf(fid,'bby = %5.4f;\n', bby);
                fprintf(fid,'               % additional capital endowment banks receive every period, as a percentage of output.\n');
                fprintf(fid,'               % (This is used for the experiments titled\n');
                fprintf(fid,'               % "Economy with More Bank Capital"): for that, set bby to  0.007\n');
                fprintf(fid,'nu = %15.14f;\n',nu);
                fprintf(fid,'% Parameters related to shock processes\n');
                fprintf(fid,'rhoz = %15.14f;    %autocorrelation of technology\n',rhoz);
                fprintf(fid,'sigmaz = %15.14f; % s.d. of innovation to technology\n',sigmaz);
                fprintf(fid,'%sigmaz = 0.0015;\n');
                fprintf(fid,'%Parameter related to capital utilisation:\n');
                fprintf(fid,'sigma_a = %15.14f;\n',sigma_a);
                fprintf(fid,'rhomp = %15.14f;  %autocorrelation of disturbances to the MP rule\n',rhomp);
                fprintf(fid,'sigmamp = %15.14f; % s.d. of innovations to these disturbances.\n',sigmamp);
                fprintf(fid,'rhobk = %15.14f;  %autocorrelation of shocks to bank capital\n',rhobk);
                fprintf(fid,'sigmabk = %15.14f; % s.d. of innovations to these disturbances.\n',rhobk);
                fprintf(fid,'%sigmabk = 0.000001;\n');
                fprintf(fid,'alpha_ss = %15.14f;\n',alpha_ss);
                fprintf(fid,'Chi =  %15.13f;            % 3.0792;       %sensitivity of entrepreneurial private benefit to monitoring intensity\n',Chi);
                fprintf(fid,'epsb = %15.13f;      %linkage between shirking and premium paid\n', epsb);
                fprintf(fid,'Btopbar =  %15.13f; %maximum private benefits from shirking;\n', Btopbar);
                fprintf(fid,'Blowbar = %5.4f; \n',Blowbar);
                fprintf(fid,'omega = %15.14f;     % Strength of response of monetary authority to deviations from credittogdp \n', omega);
                fprintf(fid,'varsigma = %15.14f;  % Strength of the endogenous link between credit to GDP devation from SS to banking sector riskiness\n', varsigma);
                fprintf(fid,'mu_ss = %15.14f;\n',mu_ss);
                fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');  
                
   fprintf(fid,'%%%%%%%%%%%%%% \n');
   fprintf(fid,'alphatemp = alpha_ss; \n');
   fprintf(fid,'mutemp = mu_ss; \n');
   fprintf(fid,'smallbtemp = Btopbar*(1+Chi*mutemp)^(-epsb); \n');
   fprintf(fid,'kk1 = delta*theta_k/(alphatemp*bigR*(1/bet-1+delta));  \n');
   fprintf(fid,'betatemp  = 1/bet; \n');
   fprintf(fid,'infl_temp = pi_ss; \n');
   fprintf(fid,'qtemp = kk1*(alphatemp*tau_b*mutemp*betatemp/delalpha-(1+mutemp+alphatemp*mutemp/(Rd_temp*delalpha))); \n');
   fprintf(fid,'qtemp = qtemp/(alphatemp*smallbtemp*kk1/(Rd_temp*delalpha)-alphatemp*bigR*kk1/Rd_temp-theta_e-theta_b-bby-alphatemp*tau_e*smallbtemp*kk1*betatemp/delalpha); \n');
   fprintf(fid,'Gtemp = 1+mutemp-(qtemp*alphatemp/Rd_temp)*(bigR-mutemp/(qtemp*delalpha)-smallbtemp/delalpha); \n');
   fprintf(fid,'IY = kk1/qtemp; \n');
   fprintf(fid,'KY = alphatemp*bigR*IY/delta;  \n');
   fprintf(fid,'KeY = tau_e*alphatemp*smallbtemp*IY/delalpha; \n');
   fprintf(fid,'KbY = tau_b*alphatemp*mutemp*IY/(qtemp*delalpha); \n');
   fprintf(fid,'KhY = KY - KeY - KbY; \n');
   fprintf(fid,'CeY = (1-tau_e)*alphatemp*smallbtemp*IY*qtemp/delalpha; \n');
   fprintf(fid,'CbY = nu*(1-tau_b)*alphatemp*mutemp*IY/delalpha; \n');
   fprintf(fid,'ChY = 1 + bby - CeY - CbY  - (1+mutemp)*IY; \n');
   fprintf(fid,'essai = (1-bet*habit)*theta_h/((1-habit)*ChY*psi_l*mark_w); \n');
   fprintf(fid,'smallh = essai^(1/(1+elas_l)); \n');
   fprintf(fid,'Htemp = smallh*eta_h^(xi_w/(xi_w-1)); \n');
   fprintf(fid,'kk2 = (KY^(1/(1-theta_k))) * (eta_e^(theta_e/(1-theta_k)))* (eta_b^(theta_b/(1-theta_k)))* mark_p^(1/(theta_k-1)); \n');
   fprintf(fid,'Ktemp = Htemp^(theta_h/(1-theta_k))* kk2; \n');
   fprintf(fid,'Ytemp = (1/mark_p)* (Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b); \n');
   fprintf(fid,'bigtheta = (mark_p-1)*Ytemp; \n');
   fprintf(fid,'Itemp = IY*Ytemp; \n');
   fprintf(fid,'rktemp = (1/mark_p)*theta_k*(Ktemp^(theta_k-1))*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^theta_b); \n');
   fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%% \n');
   fprintf(fid,'chi1 = rktemp; \n');
   fprintf(fid,'chi2 = sigma_a*chi1;   \n');
   fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%% \n');
   fprintf(fid,'w_htemp = (1/mark_p)*theta_h*(Ktemp^theta_k)*(Htemp^(theta_h-1))*(eta_e^theta_e)*(eta_b^theta_b); \n');
   fprintf(fid,'w_etemp = (1/mark_p)*theta_e*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^(theta_e-1))*(eta_b^theta_b); \n');
   fprintf(fid,'w_btemp = (1/mark_p)*theta_b*(Ktemp^theta_k)*(Htemp^theta_h)*(eta_e^theta_e)*(eta_b^(theta_b-1)); \n');
   fprintf(fid,'Ketemp = KeY*Ytemp; \n');
   fprintf(fid,'Kbtemp = KbY*Ytemp; \n');
   fprintf(fid,'Khtemp = KhY*Ytemp; \n');
   fprintf(fid,'Cetemp = CeY*Ytemp; \n');
   fprintf(fid,'Cbtemp = CbY*Ytemp; \n');
   fprintf(fid,'chtemp = ChY*Ytemp/eta_h; \n');
   fprintf(fid,'bigAtemp = eta_b*w_btemp+(rktemp+qtemp*(1-delta))*Kbtemp+bby*Ytemp; \n');
   fprintf(fid,'bigNtemp = eta_e*w_etemp+(rktemp+qtemp*(1-delta))*Ketemp; \n');
   fprintf(fid,'totCtemp = Cetemp+Cbtemp+eta_h*chtemp; \n');
   fprintf(fid,'lamtemp = (1-bet*habit)/((1-habit)*(chtemp)); \n');
   fprintf(fid,'smalldtemp = alphatemp*qtemp*( bigR - smallbtemp/delalpha - mutemp/(qtemp*delalpha))*Itemp/(Rd_temp*eta_b); \n');
   fprintf(fid,'ptemp = (Rd_temp-1)*lamtemp/ ( (Rd_temp-1)*lamtemp*eta_b*smalldtemp + eta_h*reta); \n');
   fprintf(fid,'mctemp = reta*ptemp/((Rd_temp-1)*lamtemp); \n');
   fprintf(fid,'Ratemp = alphatemp*qtemp*(mutemp/(qtemp*delalpha))*Itemp/bigAtemp; \n');
   fprintf(fid,'TLtemp = Itemp-bigNtemp; \n');
   fprintf(fid,'wtildetemp = w_htemp*infl_temp*eta_h^(-1/(1-xi_w)); \n');
   fprintf(fid,'numwtemp =  ( ( (w_htemp*infl_temp)^(xi_w*(1+elas_l)) * Htemp^(1+elas_l) )/(1-bet*phi_w) )^(1/(1+elas_l*xi_w)); \n');
   fprintf(fid,'denwtemp = ( (w_htemp^(xi_w)*infl_temp^(xi_w-1)*Htemp*lamtemp)/(1-bet*phi_w) )^(1/(1+elas_l*xi_w)); \n');
   fprintf(fid,'wtildetemp = (mark_w*psi_l)^(1/(1+xi_w*elas_l))*numwtemp/denwtemp; \n');
   fprintf(fid,'stemp = 1/mark_p; \n');
   fprintf(fid,'mgrowthtemp = infl_temp; \n');
   fprintf(fid,'expinfltemp = infl_temp; \n');
   fprintf(fid,'gYtemp = 1; \n');
   fprintf(fid,'utemp = 1; \n');
   fprintf(fid,'gammagtemp = (Itemp-bigNtemp )/bigAtemp; \n');
   fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'// 3. Model declaration \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'model;   \n');
   fprintf(fid,'1/(ch - habit*ch(-1)) - bet*habit*1/(ch(+1) - habit*ch) = lam; \n');
   fprintf(fid,'reta/(mc/p) = (Rd-1)*lam; \n');
   fprintf(fid,'rk = chi1 + chi2*(u - 1); \n');
   fprintf(fid,'q*lam = bet* (lam(+1)* ( q(+1)*(1-delta) + rk(+1)*u- ( chi1*(u-1)+0.5*chi2*(u(+1)-1)^2 ) )); \n');
   fprintf(fid,'lam = bet*lam(+1)*Rd(+1)/infl(+1); \n');
   fprintf(fid,'infl = p*mgrowth/p(-1); \n');
   fprintf(fid,'wtilde = (xi_w *psi_l/(xi_w-1))^(1/(1+elas_l*xi_w)) * numw /  denw; \n');
   fprintf(fid,'numw =  ( bet*phi_w*(numw(+1)^(1+elas_l*xi_w)) + (H^(1+elas_l)) *(w_h*infl)^(xi_w*(1+elas_l)) )^(1/(1+elas_l*xi_w)); \n');
   fprintf(fid,'denw =  ( bet*phi_w*(denw(+1)^(1+elas_l*xi_w)) + lam*H*(w_h^xi_w)*(infl^(xi_w-1)) )^(1/(1+elas_l*xi_w)) ; \n');
   fprintf(fid,'infl^(1-xi_w) *w_h^(1-xi_w) = eta_h*(1-phi_w)*wtilde^(1-xi_w) + phi_w* infl(-1)^(1-xi_w) *w_h(-1)^(1-xi_w); \n');
   fprintf(fid,'ptilde = xi_p * nump/ ((xi_p - 1) * denp ); \n');
   fprintf(fid,'nump = bet*phi_p*nump(+1) + lam*s*Y*infl^xi_p; \n');
   fprintf(fid,'denp = bet*phi_p*denp(+1) + lam*Y*infl^(xi_p-1); \n');
   fprintf(fid,'infl^(1-xi_p) = (1-phi_p)*ptilde^(1-xi_p) + phi_p*infl(-1)^(1-xi_p); \n');
   fprintf(fid,'rk = s * exp(lz)* theta_k * keff^(theta_k-1) * H^theta_h * eta_e^theta_e *eta_b^(theta_b); \n');
   fprintf(fid,'w_h = s * exp(lz)* theta_h * keff^(theta_k) * H^(theta_h-1) * eta_e^theta_e * eta_b^(theta_b); \n');
   fprintf(fid,'w_e = s * exp(lz)* theta_e * keff^(theta_k) * H^(theta_h) * eta_e^(theta_e-1) * eta_b^(theta_b); \n');
   fprintf(fid,'w_b = s * exp(lz)* theta_b * keff^(theta_k) * H^(theta_h) * eta_e^(theta_e) * eta_b^(theta_b-1); \n');
   fprintf(fid,'bigA*Ra = alpha * mu*I/delalpha; \n');
   fprintf(fid,'Rd*(1-eta_h*mc)/p = q*alpha*( bigR-smallb/delalpha - mu/(q*delalpha))*I; \n');
   fprintf(fid,'Kb = tau_b*alpha*mu*I/(q*delalpha); \n');
   fprintf(fid,'Ke = tau_e*alpha*smallb*I/delalpha; \n');
   fprintf(fid,'bigA = exp(lbk)*( rk + q*(1-delta))*Kb(-1) + eta_b*w_b; \n');
   fprintf(fid,'bigN = ( rk + q*(1-delta))*Ke(-1) + eta_e*w_e; \n');
   fprintf(fid,'Cb = (1-tau_b)*alpha*mu*I/delalpha; \n');
   fprintf(fid,'Ce = (1-tau_e)*alpha*q*smallb*I/delalpha; \n');
   fprintf(fid,'eta_h *ch + Ce + Cb + I + I*mu = Y; \n');
   fprintf(fid,'K = (1-delta)*K(-1) + alpha*bigR*I; \n');
   fprintf(fid,'log(rnot/steady_state(rnot)) = lam_r*log( rnot(-1)/steady_state(rnot) ) + (1-lam_r) * ( lam_pi*log(infl/steady_state(infl)) + lam_y*log(Y/steady_state(Y)) ) + lmp; \n');
   fprintf(fid,'Rd = rnot; \n');
   fprintf(fid,'Y = exp(lz)* (keff^theta_k)*(H^theta_h)*(eta_e^theta_e)*(eta_b^theta_b) - bigtheta; \n');
   fprintf(fid,'eta_b* smalld = (1-eta_h*mc)/p; \n');
   fprintf(fid,'keff = u*( K(-1) - Ke(-1) - Kb(-1) ) + Ke(-1) + Kb(-1); \n');
   fprintf(fid,'expinfl = infl(+1); \n');
   fprintf(fid,'gY = Y/Y(-1); \n');
   fprintf(fid,'TL = I - bigN; \n');
   fprintf(fid,'totC = eta_h *ch + Ce + Cb; \n');
   fprintf(fid,'lz = rhoz*lz(-1)+ z_shk; \n');
   fprintf(fid,'lmp = rhomp*lmp(-1)+ mp_shk; \n');
   fprintf(fid,'lbk = rhobk*lbk(-1)- bk_shk; \n');
   fprintf(fid,'log_y = log(Y); \n');
   fprintf(fid,'log_I = log(I); \n');
   fprintf(fid,'smallb = Btopbar *(1+Chi*mu)^(-epsb); \n');
   fprintf(fid,'alpha = 1 - (1-alpha_ss) - ((I-bigN)/Y-(steady_state(I)-steady_state(bigN))/steady_state(Y))*varsigma ;  \n');
   fprintf(fid,'gammag = steady_state(gammag) + omega*((I-bigN-CBBS*Y)/Y - (steady_state(I)-steady_state(bigN))/steady_state(Y));  \n');
   fprintf(fid,'I =gammag * bigA + bigN + CBBS*Y; \n');
   fprintf(fid,'gammag = 1 + (q*Ra/Rd)*(delalpha*bigR-smallb-mu/q)/mu - Ra*delalpha/alpha;  \n');
   fprintf(fid,'CBBS = lqe/steady_state(Y); \n');
   fprintf(fid,'lqe = rhoqe*lqe(-1)+ qe_shk; \n');
   fprintf(fid,'bdp = mu/q; \n');
   fprintf(fid,'inflexp1  = infl(+1); inflexp2  = inflexp1(+1); inflexp3  = inflexp2(+1); inflexp4  = inflexp3(+1); inflexp5  = inflexp4(+1); inflexp6  = inflexp5(+1); inflexp7  = inflexp6(+1); inflexp8  = inflexp7(+1); inflexp9  = inflexp8(+1); inflexp10  = inflexp9(+1); inflexp11  = inflexp10(+1); inflexp12  = inflexp11(+1); inflexp13  = inflexp12(+1); inflexp14  = inflexp13(+1); inflexp15  = inflexp14(+1); inflexp16  = inflexp15(+1); inflexp17  = inflexp16(+1); inflexp18  = inflexp17(+1); inflexp19  = inflexp18(+1); inflexp20  = inflexp19(+1); inflexp21  = inflexp20(+1); inflexp22  = inflexp21(+1); inflexp23  = inflexp22(+1); inflexp24  = inflexp23(+1); inflexp25  = inflexp24(+1); inflexp26  = inflexp25(+1); inflexp27  = inflexp26(+1); inflexp28  = inflexp27(+1); inflexp29  = inflexp28(+1); inflexp30  = inflexp29(+1); inflexp31  = inflexp30(+1); inflexp32  = inflexp31(+1); inflexp33  = inflexp32(+1); inflexp34  = inflexp33(+1); inflexp35  = inflexp34(+1); inflexp36  = inflexp35(+1); inflexp37  = inflexp36(+1); inflexp38  = inflexp37(+1); inflexp39  = inflexp38(+1); inflexp40  = inflexp39(+1); Rexp1  = Rd(+1); Rexp2  = Rexp1(+1); Rexp3  = Rexp2(+1); Rexp4  = Rexp3(+1); Rexp5  = Rexp4(+1); Rexp6  = Rexp5(+1); Rexp7  = Rexp6(+1); Rexp8  = Rexp7(+1); Rexp9  = Rexp8(+1); Rexp10  = Rexp9(+1); Rexp11  = Rexp10(+1); Rexp12  = Rexp11(+1); Rexp13  = Rexp12(+1); Rexp14  = Rexp13(+1); Rexp15  = Rexp14(+1); Rexp16  = Rexp15(+1); Rexp17  = Rexp16(+1); Rexp18  = Rexp17(+1); Rexp19  = Rexp18(+1); Rexp20  = Rexp19(+1); Rexp21  = Rexp20(+1); Rexp22  = Rexp21(+1); Rexp23  = Rexp22(+1); Rexp24  = Rexp23(+1); Rexp25  = Rexp24(+1); Rexp26  = Rexp25(+1); Rexp27  = Rexp26(+1); Rexp28  = Rexp27(+1); Rexp29  = Rexp28(+1); Rexp30  = Rexp29(+1); Rexp31  = Rexp30(+1); Rexp32  = Rexp31(+1); Rexp33  = Rexp32(+1); Rexp34  = Rexp33(+1); Rexp35  = Rexp34(+1); Rexp36  = Rexp35(+1); Rexp37  = Rexp36(+1); Rexp38  = Rexp37(+1); Rexp39  = Rexp38(+1); Rexp40  = Rexp39(+1); \n');
   fprintf(fid,'util = log( ch - habit*ch(-1) ) -psi_l * (1-bet*habit)*theta_h/((1-habit)*ch*psi_l*mark_w)  / (1+elas_l) + reta*log(mc/p); \n');
   fprintf(fid,'welf = util + bet*util(+1); \n');
   fprintf(fid,'end; \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'// 4. Initial values and steady state \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'resid; \n');
   fprintf(fid,'steady;   \n');
   fprintf(fid,'check;                     \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'// 5. Shock declaration   \n');
   fprintf(fid,'//                        \n');
   fprintf(fid,'//--------------------------------------------------------------------- \n');
   fprintf(fid,'shocks; \n');
   fprintf(fid,'var z_shk = sigmaz^2; \n');
   fprintf(fid,'var mp_shk = sigmamp^2; \n');
   fprintf(fid,'var bk_shk = sigmabk^2; \n');
   fprintf(fid,'var qe_shk = sigmaqe^2; \n');
   fprintf(fid,'end; \n');
   fprintf(fid,'%stoch_simul(order=2,irf = 40,hp_filter = 1600)  Y I q gammag bigA TL bigN Rd infl mu alpha; \n');
   fprintf(fid,'%options_.debug = 1; \n');
   fprintf(fid,'optim_weights; \n');
   fprintf(fid,'welf 1; \n');
   fprintf(fid,'end; \n');
   fprintf(fid,'osr_params omega; \n');
   fprintf(fid,'%osr_params_bounds; \n');
   fprintf(fid,'%omega, -100 , 0; \n');
   fprintf(fid,'%end; \n');
   fprintf(fid,'osr(opt_algo = 9, nograph); \n');