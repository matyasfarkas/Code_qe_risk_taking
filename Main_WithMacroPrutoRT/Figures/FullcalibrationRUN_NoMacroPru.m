%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

clear all
close all

addpath(('C:\dynare\4.5.7\matlab'));
%% Parameters to change
elas = -1.27951709303868;  % Estimated elasticity from the NSR identification
maxiter = 100;
lam_rr = zeros(maxiter,1);
lam_rr(1)= 0.8;        % AR(1) parameter in the Taylor rule, for short term interest rate smoothing

lam_pir = zeros(maxiter,1);
lam_pir(1) = 1.5;       % Parameter on the inflation target in the Taylor rule

lam_yr = zeros(maxiter,1);
lam_yr(1) = 0;      % Parameter of output stabilization in the Taylor rule

lambdaa = 2;       % Newton Rhapson root search parameter
est_elas = zeros(maxiter,1);

mkdir('RP');
copyfile RP_omega_write.m RP
copyfile Parameters_calibration.m RP

mkdir('Calibration_tmp');
copyfile writefiles_nomacropru.m Calibration_tmp
copyfile Parameters_calibration.m Calibration_tmp
copyfile setpathdynare442.m Calibration_tmp

lam_rr = [0.9];
lam_yr = [0.0:0.05:0.1];
lam_pir = [1.7:0.05:1.9];
CBBSpr = [-0.5:-0.5:-1.5];
for jjjj = 1:3
for jjj = 1:4 % 1:2
    
    for jj = 1 %1:2
        
        for j =  1 %1:2
            %%
            Parameters_calibration
            CBBSp = CBBSpr(jjjj);
            if j == 1 && jj ==1 && jjj==1
                CBBSp = CBBSpr(jjjj);
                lam_r = lam_rr(j);
                lam_y = lam_yr(jj);
                lam_pi = lam_pir(jjj);
                
                save rootfinder.mat lam_r  lam_pi lam_y lambdaa maxiter elas est_elas CBBSp
            else
                load rootfinder.mat
                %
                %             j = j +1;
                %             jj = jj+1;
                %             jjj= jjj+1;
            end
            lam_r = lam_rr(j);
            lam_y = lam_yr(jj);
            lam_pi = lam_pir(jjj);
            CBBSp = CBBSpr(jjjj);

            %%
            try
                cd RP
                RP_omega_write
                addpath(('C:\dynare\4.5.7\matlab'));
                dynare RP_omega.mod noclearall;
                omega = oo_.osr.optim_params.omega;
                
                cd ..
                cd Calibration_tmp
                setpathdynare442
                writefiles_nomacropru
                
                global M_ oo_
                
                % modnam and modnamstar below choose model
                modnam = 'FMcc'; % non binding baseline model
                modnamstar = 'FMzlb'; % Model with ZLB and QE
                modnamstar_noqe = 'noqeFMzlb'; % Model with ZLB and wo QE
                modnamstar_nomacropru = 'nomacropruFMzlb'; % Model with ZLB QE and wo MacroPru
                modnam_nomacropru  = 'FMtiv'; % non binding baseline model
                
                % see notes 1 and 2 in readme file
                constraint = 'Rd< -(pi_ss/bet)+1';
                constraint_relax ='rnot>-(pi_ss/bet)+1';
                
                
                % Pick innovation for IRFs
                irfshock =char('bk_shk');      % label for innovation for IRFs
                % needs to be an exogenous variable in the
                % dynare .mod files
                
                
                % Shock path
                nper=10;
                
                shockssequence = [
                    1*ones(nper,1)/nper*2
                    ];         % scale factor for simulations
                nperiods = size(shockssequence,1)+30;            %length of IRFs
                
                maxiter = 10000;
                
                %% Solve model
                
                % Solve QE model against baseline, generate model IRFs
                tic
                [zdatalinear zdatapiecewise zdatass oobase_ Mbase_  ] = ...
                    solve_one_constraint(modnam,modnamstar,...
                    constraint, constraint_relax,...
                    shockssequence,irfshock,nperiods,maxiter);
                toc
                
                %%             % Solve no macropru model against baseline, generate model IRFs
                constraint2 = 'gammag >2.5';
                constraint2_relax ='gammag < 2.5';
                
                [zdatalinear zdatapiecewise_nomacropru zdatass2 oo2 M_2] = ...
                    solve_two_constraints(modnam,modnamstar,modnamstar_nomacropru,modnamstar_nomacropru,...
                    constraint, constraint2,constraint_relax, constraint2_relax,...
                    shockssequence,irfshock,nperiods,1,maxiter);
                %%
                % Solve no macropru model against baseline, generate model IRFs
%                 tic
%                 [zdatapiecewise_tiv zdatapiecewise_nomacropru zdatass1 oobase1_ Mbase1_  ] = ...
%                     solve_one_constraint(modnam,modnamstar_nomacropru,...
%                     constraint, constraint_relax,...
%                     shockssequence,irfshock,nperiods,maxiter);
%                 toc
                
                
                % Solve no QE no ZLB model against baseline, generate model IRFs
                % Solve model, generate model IRFs
                tic
                [zdatalinear zdatapiecewise_noqe zdatass1 oobase1_ Mbase1_  ] = ...
                    solve_one_constraint(modnam,modnamstar_noqe,...
                    constraint, constraint_relax,...
                    shockssequence,irfshock,nperiods,maxiter);
                toc
                
                % unpack the IRFs
                for i=1:M_.endo_nbr
                    eval([deblank(M_.endo_names(i,:)),'_uncdifference=zdatalinear(:,i);']);
                    eval([deblank(M_.endo_names(i,:)),'_difference=zdatapiecewise(:,i);']);
                    eval([deblank(M_.endo_names(i,:)),'_difference_nomacropru=zdatapiecewise_nomacropru(:,i);']);
                    eval([deblank(M_.endo_names(i,:)),'_difference_noqe=zdatapiecewise_noqe(:,i);']);
                    eval([deblank(M_.endo_names(i,:)),'_ss=zdatass(i);']);
                end
                
                
                % get parameter values
                for i=1:Mbase_.param_nbr
                    eval([Mbase_.param_names(i,:),'=Mbase_.params(i);'])
                end
                
                %% Modify to plot IRFs
                
                titlelist = char('Rd, Policy Interest Rate (percent at AR)','bigA, Bank Capital',...
                    'I, Investment','infl, Inflatiom', 'CBBS, QE', ...
                    'q, Real price of capital', ...
                    'VaR, Bank default probability', ...
                    'gammag, Leverage ratio', ...
                    'Rb/R, Share of banker''s profits from overall', ...
                    'expinfl, Inflation expectations');
                
                percent = '  % dev. from SS    ';
                level =   '      Level         ';
                ylabels = char(level,... % Policy rate
                   percent,... % Bank capital rate
                   percent,... % Investments
                   level,...  % inflation
                   percent,... % QE
                   percent,... % q
                   percent,...  % Var
                   level,... % Gammag
                   percent,... % Rb/R
                   level); % inflation expectations
                
                figtitle = '';
                % ZLB with QE
                line1=[400*(Rd_difference+Rd_ss-1),bigA_difference,I_difference,400*(infl_difference+pi_ss-1),CBBS_difference,q_difference,-BVaR_difference,1./(8+gammag_difference),Ra_difference/bigR*100,(expinfl_difference+pi_ss-1)*400];
                % No ZLB no QE
                line2=[400*(Rd_uncdifference+Rd_ss-1),bigA_uncdifference,I_uncdifference,400*(infl_uncdifference+pi_ss-1),CBBS_uncdifference,q_uncdifference,-BVaR_uncdifference,1./(8+gammag_uncdifference),Ra_uncdifference/bigR*100,(expinfl_uncdifference+pi_ss-1)*400];
                %ZLB, QE, no MacroPru
                line3=[400*(Rd_difference_nomacropru+Rd_ss-1),bigA_difference_nomacropru,I_difference_nomacropru,400*(infl_difference_nomacropru+pi_ss-1),CBBS_difference_nomacropru,q_difference_nomacropru,-BVaR_difference_nomacropru,1./(8+gammag_difference_nomacropru),Ra_difference_nomacropru/bigR*100,(expinfl_difference_nomacropru+pi_ss-1)*400];
                %ZLB, no QE with Macropu
                line4=[400*(Rd_difference_noqe+Rd_ss-1),bigA_difference_noqe,I_difference_noqe,400*(infl_difference_noqe+pi_ss-1),CBBS_difference_noqe,q_difference_noqe,-BVaR_difference_noqe,1./(8+gammag_difference_noqe),Ra_difference_noqe/bigR*100,(expinfl_difference_noqe+pi_ss-1)*400];
                
                
                
                legendlist = cellstr(char('ZLB with QE','ZLB, no QE','ZLB, TIV MacroPru'));
                figlabel = '';
                makechart(titlelist,legendlist,figlabel,ylabels,line1,line4,line3)
                diff_infexp = line1(:,end)-line4(:,end);
                diff_rates = line1(:,1)-line4(:,1);
                
                roundn = @(x,n) round(x.*10.^n)./10.^n;
                for ii = 1:40
                    eval(['inflexp', num2str(ii), '_difference = roundn(inflexp', num2str(ii), '_difference ,10);']);
                    eval(['inflexp', num2str(ii), '_uncdifference = roundn(inflexp', num2str(ii), '_uncdifference ,10);']);
                    eval(['inflexp', num2str(ii), '_difference_nomacropru = roundn(inflexp', num2str(ii), '_difference_nomacropru ,10);']);
                    eval(['inflexp', num2str(ii), '_difference_noqe = roundn(inflexp', num2str(ii), '_difference_noqe ,10);']);
                    
                end
                
                
                inf40_difference = (((pi_ss+inflexp1_difference).*(pi_ss+inflexp2_difference).*(pi_ss+inflexp3_difference).*(pi_ss+inflexp4_difference).*(pi_ss+inflexp5_difference).*(pi_ss+inflexp6_difference).*(pi_ss+inflexp7_difference).*(pi_ss+inflexp8_difference).*(pi_ss+inflexp9_difference).*(pi_ss+inflexp10_difference).*(pi_ss+inflexp11_difference).*(pi_ss+inflexp12_difference).*(pi_ss+inflexp13_difference).*(pi_ss+inflexp14_difference).*(pi_ss+inflexp15_difference).*(pi_ss+inflexp16_difference).*(pi_ss+inflexp17_difference).*(pi_ss+inflexp18_difference).*(pi_ss+inflexp19_difference).*(pi_ss+inflexp20_difference).*(pi_ss+inflexp21_difference).*(pi_ss+inflexp22_difference).*(pi_ss+inflexp23_difference).*(pi_ss+inflexp24_difference).*(pi_ss+inflexp25_difference).*(pi_ss+inflexp26_difference).*(pi_ss+inflexp27_difference).*(pi_ss+inflexp28_difference).*(pi_ss+inflexp29_difference).*(pi_ss+inflexp30_difference).*(pi_ss+inflexp31_difference).*(pi_ss+inflexp32_difference).*(pi_ss+inflexp33_difference).*(pi_ss+inflexp34_difference).*(pi_ss+inflexp35_difference).*(pi_ss+inflexp36_difference).*(pi_ss+inflexp37_difference).*(pi_ss+inflexp38_difference).*(pi_ss+inflexp39_difference).*(pi_ss+inflexp40_difference)).^(1/40)-1)*400;
                inf40_uncdifference= (((pi_ss+inflexp1_uncdifference).*(pi_ss+inflexp2_uncdifference).*(pi_ss+inflexp3_uncdifference).*(pi_ss+inflexp4_uncdifference).*(pi_ss+inflexp5_uncdifference).*(pi_ss+inflexp6_uncdifference).*(pi_ss+inflexp7_uncdifference).*(pi_ss+inflexp8_uncdifference).*(pi_ss+inflexp9_uncdifference).*(pi_ss+inflexp10_uncdifference).*(pi_ss+inflexp11_uncdifference).*(pi_ss+inflexp12_uncdifference).*(pi_ss+inflexp13_uncdifference).*(pi_ss+inflexp14_uncdifference).*(pi_ss+inflexp15_uncdifference).*(pi_ss+inflexp16_uncdifference).*(pi_ss+inflexp17_uncdifference).*(pi_ss+inflexp18_uncdifference).*(pi_ss+inflexp19_uncdifference).*(pi_ss+inflexp20_uncdifference).*(pi_ss+inflexp21_uncdifference).*(pi_ss+inflexp22_uncdifference).*(pi_ss+inflexp23_uncdifference).*(pi_ss+inflexp24_uncdifference).*(pi_ss+inflexp25_uncdifference).*(pi_ss+inflexp26_uncdifference).*(pi_ss+inflexp27_uncdifference).*(pi_ss+inflexp28_uncdifference).*(pi_ss+inflexp29_uncdifference).*(pi_ss+inflexp30_uncdifference).*(pi_ss+inflexp31_uncdifference).*(pi_ss+inflexp32_uncdifference).*(pi_ss+inflexp33_uncdifference).*(pi_ss+inflexp34_uncdifference).*(pi_ss+inflexp35_uncdifference).*(pi_ss+inflexp36_uncdifference).*(pi_ss+inflexp37_uncdifference).*(pi_ss+inflexp38_uncdifference).*(pi_ss+inflexp39_uncdifference).*(pi_ss+inflexp40_uncdifference)).^(1/40)-1)*400;
                inf40_difference_nomacropru= (((pi_ss+inflexp1_difference_nomacropru).*(pi_ss+inflexp2_difference_nomacropru).*(pi_ss+inflexp3_difference_nomacropru).*(pi_ss+inflexp4_difference_nomacropru).*(pi_ss+inflexp5_difference_nomacropru).*(pi_ss+inflexp6_difference_nomacropru).*(pi_ss+inflexp7_difference_nomacropru).*(pi_ss+inflexp8_difference_nomacropru).*(pi_ss+inflexp9_difference_nomacropru).*(pi_ss+inflexp10_difference_nomacropru).*(pi_ss+inflexp11_difference_nomacropru).*(pi_ss+inflexp12_difference_nomacropru).*(pi_ss+inflexp13_difference_nomacropru).*(pi_ss+inflexp14_difference_nomacropru).*(pi_ss+inflexp15_difference_nomacropru).*(pi_ss+inflexp16_difference_nomacropru).*(pi_ss+inflexp17_difference_nomacropru).*(pi_ss+inflexp18_difference_nomacropru).*(pi_ss+inflexp19_difference_nomacropru).*(pi_ss+inflexp20_difference_nomacropru).*(pi_ss+inflexp21_difference_nomacropru).*(pi_ss+inflexp22_difference_nomacropru).*(pi_ss+inflexp23_difference_nomacropru).*(pi_ss+inflexp24_difference_nomacropru).*(pi_ss+inflexp25_difference_nomacropru).*(pi_ss+inflexp26_difference_nomacropru).*(pi_ss+inflexp27_difference_nomacropru).*(pi_ss+inflexp28_difference_nomacropru).*(pi_ss+inflexp29_difference_nomacropru).*(pi_ss+inflexp30_difference_nomacropru).*(pi_ss+inflexp31_difference_nomacropru).*(pi_ss+inflexp32_difference_nomacropru).*(pi_ss+inflexp33_difference_nomacropru).*(pi_ss+inflexp34_difference_nomacropru).*(pi_ss+inflexp35_difference_nomacropru).*(pi_ss+inflexp36_difference_nomacropru).*(pi_ss+inflexp37_difference_nomacropru).*(pi_ss+inflexp38_difference_nomacropru).*(pi_ss+inflexp39_difference_nomacropru).*(pi_ss+inflexp40_difference_nomacropru)).^(1/40)-1)*400;
                inf40_difference_noqe= (((pi_ss+inflexp1_difference_noqe).*(pi_ss+inflexp2_difference_noqe).*(pi_ss+inflexp3_difference_noqe).*(pi_ss+inflexp4_difference_noqe).*(pi_ss+inflexp5_difference_noqe).*(pi_ss+inflexp6_difference_noqe).*(pi_ss+inflexp7_difference_noqe).*(pi_ss+inflexp8_difference_noqe).*(pi_ss+inflexp9_difference_noqe).*(pi_ss+inflexp10_difference_noqe).*(pi_ss+inflexp11_difference_noqe).*(pi_ss+inflexp12_difference_noqe).*(pi_ss+inflexp13_difference_noqe).*(pi_ss+inflexp14_difference_noqe).*(pi_ss+inflexp15_difference_noqe).*(pi_ss+inflexp16_difference_noqe).*(pi_ss+inflexp17_difference_noqe).*(pi_ss+inflexp18_difference_noqe).*(pi_ss+inflexp19_difference_noqe).*(pi_ss+inflexp20_difference_noqe).*(pi_ss+inflexp21_difference_noqe).*(pi_ss+inflexp22_difference_noqe).*(pi_ss+inflexp23_difference_noqe).*(pi_ss+inflexp24_difference_noqe).*(pi_ss+inflexp25_difference_noqe).*(pi_ss+inflexp26_difference_noqe).*(pi_ss+inflexp27_difference_noqe).*(pi_ss+inflexp28_difference_noqe).*(pi_ss+inflexp29_difference_noqe).*(pi_ss+inflexp30_difference_noqe).*(pi_ss+inflexp31_difference_noqe).*(pi_ss+inflexp32_difference_noqe).*(pi_ss+inflexp33_difference_noqe).*(pi_ss+inflexp34_difference_noqe).*(pi_ss+inflexp35_difference_noqe).*(pi_ss+inflexp36_difference_noqe).*(pi_ss+inflexp37_difference_noqe).*(pi_ss+inflexp38_difference_noqe).*(pi_ss+inflexp39_difference_noqe).*(pi_ss+inflexp40_difference_noqe)).^(1/40)-1)*400;
                
                
                R40_difference = (((Rd_ss+Rexp1_difference).*(Rd_ss+Rexp2_difference).*(Rd_ss+Rexp3_difference).*(Rd_ss+Rexp4_difference).*(Rd_ss+Rexp5_difference).*(Rd_ss+Rexp6_difference).*(Rd_ss+Rexp7_difference).*(Rd_ss+Rexp8_difference).*(Rd_ss+Rexp9_difference).*(Rd_ss+Rexp10_difference).*(Rd_ss+Rexp11_difference).*(Rd_ss+Rexp12_difference).*(Rd_ss+Rexp13_difference).*(Rd_ss+Rexp14_difference).*(Rd_ss+Rexp15_difference).*(Rd_ss+Rexp16_difference).*(Rd_ss+Rexp17_difference).*(Rd_ss+Rexp18_difference).*(Rd_ss+Rexp19_difference).*(Rd_ss+Rexp20_difference).*(Rd_ss+Rexp21_difference).*(Rd_ss+Rexp22_difference).*(Rd_ss+Rexp23_difference).*(Rd_ss+Rexp24_difference).*(Rd_ss+Rexp25_difference).*(Rd_ss+Rexp26_difference).*(Rd_ss+Rexp27_difference).*(Rd_ss+Rexp28_difference).*(Rd_ss+Rexp29_difference).*(Rd_ss+Rexp30_difference).*(Rd_ss+Rexp31_difference).*(Rd_ss+Rexp32_difference).*(Rd_ss+Rexp33_difference).*(Rd_ss+Rexp34_difference).*(Rd_ss+Rexp35_difference).*(Rd_ss+Rexp36_difference).*(Rd_ss+Rexp37_difference).*(Rd_ss+Rexp38_difference).*(Rd_ss+Rexp39_difference).*(Rd_ss+Rexp40_difference)).^(1/40)-1)*400;
                R40_uncdifference= (((Rd_ss+Rexp1_uncdifference).*(Rd_ss+Rexp2_uncdifference).*(Rd_ss+Rexp3_uncdifference).*(Rd_ss+Rexp4_uncdifference).*(Rd_ss+Rexp5_uncdifference).*(Rd_ss+Rexp6_uncdifference).*(Rd_ss+Rexp7_uncdifference).*(Rd_ss+Rexp8_uncdifference).*(Rd_ss+Rexp9_uncdifference).*(Rd_ss+Rexp10_uncdifference).*(Rd_ss+Rexp11_uncdifference).*(Rd_ss+Rexp12_uncdifference).*(Rd_ss+Rexp13_uncdifference).*(Rd_ss+Rexp14_uncdifference).*(Rd_ss+Rexp15_uncdifference).*(Rd_ss+Rexp16_uncdifference).*(Rd_ss+Rexp17_uncdifference).*(Rd_ss+Rexp18_uncdifference).*(Rd_ss+Rexp19_uncdifference).*(Rd_ss+Rexp20_uncdifference).*(Rd_ss+Rexp21_uncdifference).*(Rd_ss+Rexp22_uncdifference).*(Rd_ss+Rexp23_uncdifference).*(Rd_ss+Rexp24_uncdifference).*(Rd_ss+Rexp25_uncdifference).*(Rd_ss+Rexp26_uncdifference).*(Rd_ss+Rexp27_uncdifference).*(Rd_ss+Rexp28_uncdifference).*(Rd_ss+Rexp29_uncdifference).*(Rd_ss+Rexp30_uncdifference).*(Rd_ss+Rexp31_uncdifference).*(Rd_ss+Rexp32_uncdifference).*(Rd_ss+Rexp33_uncdifference).*(Rd_ss+Rexp34_uncdifference).*(Rd_ss+Rexp35_uncdifference).*(Rd_ss+Rexp36_uncdifference).*(Rd_ss+Rexp37_uncdifference).*(Rd_ss+Rexp38_uncdifference).*(Rd_ss+Rexp39_uncdifference).*(Rd_ss+Rexp40_uncdifference)).^(1/40)-1)*400;
                R40_difference_nomacropru= (((Rd_ss+Rexp1_difference_nomacropru).*(Rd_ss+Rexp2_difference_nomacropru).*(Rd_ss+Rexp3_difference_nomacropru).*(Rd_ss+Rexp4_difference_nomacropru).*(Rd_ss+Rexp5_difference_nomacropru).*(Rd_ss+Rexp6_difference_nomacropru).*(Rd_ss+Rexp7_difference_nomacropru).*(Rd_ss+Rexp8_difference_nomacropru).*(Rd_ss+Rexp9_difference_nomacropru).*(Rd_ss+Rexp10_difference_nomacropru).*(Rd_ss+Rexp11_difference_nomacropru).*(Rd_ss+Rexp12_difference_nomacropru).*(Rd_ss+Rexp13_difference_nomacropru).*(Rd_ss+Rexp14_difference_nomacropru).*(Rd_ss+Rexp15_difference_nomacropru).*(Rd_ss+Rexp16_difference_nomacropru).*(Rd_ss+Rexp17_difference_nomacropru).*(Rd_ss+Rexp18_difference_nomacropru).*(Rd_ss+Rexp19_difference_nomacropru).*(Rd_ss+Rexp20_difference_nomacropru).*(Rd_ss+Rexp21_difference_nomacropru).*(Rd_ss+Rexp22_difference_nomacropru).*(Rd_ss+Rexp23_difference_nomacropru).*(Rd_ss+Rexp24_difference_nomacropru).*(Rd_ss+Rexp25_difference_nomacropru).*(Rd_ss+Rexp26_difference_nomacropru).*(Rd_ss+Rexp27_difference_nomacropru).*(Rd_ss+Rexp28_difference_nomacropru).*(Rd_ss+Rexp29_difference_nomacropru).*(Rd_ss+Rexp30_difference_nomacropru).*(Rd_ss+Rexp31_difference_nomacropru).*(Rd_ss+Rexp32_difference_nomacropru).*(Rd_ss+Rexp33_difference_nomacropru).*(Rd_ss+Rexp34_difference_nomacropru).*(Rd_ss+Rexp35_difference_nomacropru).*(Rd_ss+Rexp36_difference_nomacropru).*(Rd_ss+Rexp37_difference_nomacropru).*(Rd_ss+Rexp38_difference_nomacropru).*(Rd_ss+Rexp39_difference_nomacropru).*(Rd_ss+Rexp40_difference_nomacropru)).^(1/40)-1)*400;
                R40_difference_noqe= (((Rd_ss+Rexp1_difference_noqe).*(Rd_ss+Rexp2_difference_noqe).*(Rd_ss+Rexp3_difference_noqe).*(Rd_ss+Rexp4_difference_noqe).*(Rd_ss+Rexp5_difference_noqe).*(Rd_ss+Rexp6_difference_noqe).*(Rd_ss+Rexp7_difference_noqe).*(Rd_ss+Rexp8_difference_noqe).*(Rd_ss+Rexp9_difference_noqe).*(Rd_ss+Rexp10_difference_noqe).*(Rd_ss+Rexp11_difference_noqe).*(Rd_ss+Rexp12_difference_noqe).*(Rd_ss+Rexp13_difference_noqe).*(Rd_ss+Rexp14_difference_noqe).*(Rd_ss+Rexp15_difference_noqe).*(Rd_ss+Rexp16_difference_noqe).*(Rd_ss+Rexp17_difference_noqe).*(Rd_ss+Rexp18_difference_noqe).*(Rd_ss+Rexp19_difference_noqe).*(Rd_ss+Rexp20_difference_noqe).*(Rd_ss+Rexp21_difference_noqe).*(Rd_ss+Rexp22_difference_noqe).*(Rd_ss+Rexp23_difference_noqe).*(Rd_ss+Rexp24_difference_noqe).*(Rd_ss+Rexp25_difference_noqe).*(Rd_ss+Rexp26_difference_noqe).*(Rd_ss+Rexp27_difference_noqe).*(Rd_ss+Rexp28_difference_noqe).*(Rd_ss+Rexp29_difference_noqe).*(Rd_ss+Rexp30_difference_noqe).*(Rd_ss+Rexp31_difference_noqe).*(Rd_ss+Rexp32_difference_noqe).*(Rd_ss+Rexp33_difference_noqe).*(Rd_ss+Rexp34_difference_noqe).*(Rd_ss+Rexp35_difference_noqe).*(Rd_ss+Rexp36_difference_noqe).*(Rd_ss+Rexp37_difference_noqe).*(Rd_ss+Rexp38_difference_noqe).*(Rd_ss+Rexp39_difference_noqe).*(Rd_ss+Rexp40_difference_noqe)).^(1/40)-1)*400;
                
                
                
                dR40 = R40_difference - R40_difference_noqe;
                dinf40 = inf40_difference -inf40_difference_noqe;
                index = nan;
                % searching for periods when the QE has been extending.
                for looper= 2:size(CBBS_difference,1)
                    if CBBS_difference(looper)-CBBS_difference(looper-1)>10^-3
                        if isnan(index)
                            index = looper;
                        else
                            index = [ index ; looper];
                        end
                    end
                end
                
                frame_h = get(handle(gcf),'JavaFrame');
                set(frame_h,'Maximized',1);
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                % Get rid of tool bar and pulldown menus that are along top of figure.
                set(gcf, 'Toolbar', 'none', 'Menu', 'none');
                mkdir('Figures');
                cd Figures
                saveas(gcf,['FigureFORrhor',num2str(j),'yr', num2str(jj),'pir', num2str(jjj),'CBBSp',num2str(jjjj),'.fig']);
                print(['FigureFORrhor',num2str(j),'yr', num2str(jj),'pir', num2str(jjj),'CBBSp',num2str(jjjj)],'-dpng');
                matlab2tikz(['FigureFORrhor',num2str(j),'yr', num2str(jj),'pir', num2str(jjj),'CBBSp',num2str(jjjj),'.tex']);
                est_elas = mean(dinf40(index))/mean(dR40(index));
                fid = fopen(['FigureFORrhor',num2str(j),'yr', num2str(jj),'pir', num2str(jjj),'CBBSp',num2str(jjjj),'.txt'],'w+');
                fprintf(fid,'Estimated elasticity between the long term yield and inflation expectations: \n');
                fprintf(fid,'est_elas = %5.4f;\n', est_elas);
                fclose('all');
                cd ..
                cd ..
                
            catch
                est_elas(j,jj,jjj,jjjjj) = mean(dinf40(index))/mean(dR40(index));
                mean(dinf40(index))/mean(dR40(index)),
                pause(10)
                
                
                
                if abs(est_elas(j)-elas) < 0.01
                end
                save rootfinder.mat lam_rr  lam_pir lam_yr lam_r  lam_pi lam_y lambdaa maxiter elas est_elas
                
            end
            close all
        end
        
    end
end
end
figure
plot((1:size(dR40,1)), dR40,'r--')
hold on
plot((1:size(dinf40,1)), dinf40,'g-*')
legend('10 year yields', '10 year inflation expectations')

qe = [inflexp1_difference,inflexp2_difference,inflexp3_difference,inflexp4_difference,inflexp5_difference,inflexp6_difference,inflexp7_difference,inflexp8_difference,inflexp9_difference,inflexp10_difference,inflexp11_difference,inflexp12_difference,inflexp13_difference,inflexp14_difference,inflexp15_difference,inflexp16_difference,inflexp17_difference,inflexp18_difference,inflexp19_difference,inflexp20_difference,inflexp21_difference,inflexp22_difference,inflexp23_difference,inflexp24_difference,inflexp25_difference,inflexp26_difference,inflexp27_difference,inflexp28_difference,inflexp29_difference,inflexp30_difference,inflexp31_difference,inflexp32_difference,inflexp33_difference,inflexp34_difference,inflexp35_difference,inflexp36_difference,inflexp37_difference,inflexp38_difference,inflexp39_difference,inflexp40_difference];
qe = (qe+pi_ss-1)*400;
nomacropru = [inflexp1_difference_nomacropru,inflexp2_difference_nomacropru,inflexp3_difference_nomacropru,inflexp4_difference_nomacropru,inflexp5_difference_nomacropru,inflexp6_difference_nomacropru,inflexp7_difference_nomacropru,inflexp8_difference_nomacropru,inflexp9_difference_nomacropru,inflexp10_difference_nomacropru,inflexp11_difference_nomacropru,inflexp12_difference_nomacropru,inflexp13_difference_nomacropru,inflexp14_difference_nomacropru,inflexp15_difference_nomacropru,inflexp16_difference_nomacropru,inflexp17_difference_nomacropru,inflexp18_difference_nomacropru,inflexp19_difference_nomacropru,inflexp20_difference_nomacropru,inflexp21_difference_nomacropru,inflexp22_difference_nomacropru,inflexp23_difference_nomacropru,inflexp24_difference_nomacropru,inflexp25_difference_nomacropru,inflexp26_difference_nomacropru,inflexp27_difference_nomacropru,inflexp28_difference_nomacropru,inflexp29_difference_nomacropru,inflexp30_difference_nomacropru,inflexp31_difference_nomacropru,inflexp32_difference_nomacropru,inflexp33_difference_nomacropru,inflexp34_difference_nomacropru,inflexp35_difference_nomacropru,inflexp36_difference_nomacropru,inflexp37_difference_nomacropru,inflexp38_difference_nomacropru,inflexp39_difference_nomacropru,inflexp40_difference_nomacropru];
nomacropru = (nomacropru+pi_ss-1)*400;
nozlb = [inflexp1_uncdifference,inflexp2_uncdifference,inflexp3_uncdifference,inflexp4_uncdifference,inflexp5_uncdifference,inflexp6_uncdifference,inflexp7_uncdifference,inflexp8_uncdifference,inflexp9_uncdifference,inflexp10_uncdifference,inflexp11_uncdifference,inflexp12_uncdifference,inflexp13_uncdifference,inflexp14_uncdifference,inflexp15_uncdifference,inflexp16_uncdifference,inflexp17_uncdifference,inflexp18_uncdifference,inflexp19_uncdifference,inflexp20_uncdifference,inflexp21_uncdifference,inflexp22_uncdifference,inflexp23_uncdifference,inflexp24_uncdifference,inflexp25_uncdifference,inflexp26_uncdifference,inflexp27_uncdifference,inflexp28_uncdifference,inflexp29_uncdifference,inflexp30_uncdifference,inflexp31_uncdifference,inflexp32_uncdifference,inflexp33_uncdifference,inflexp34_uncdifference,inflexp35_uncdifference,inflexp36_uncdifference,inflexp37_uncdifference,inflexp38_uncdifference,inflexp39_uncdifference,inflexp40_uncdifference];
nozlb= (nozlb+pi_ss-1)*400;

figure
mesh(nomacropru(1:50,:))
hold on
surf(nozlb(1:50,:))

ddinf = (line1(2:end,end)-line1(1:end-1,end)) -(line3(2:end,end)-line3(1:end-1,end));
ddrat = (line1(2:end,1)-line1(1:end-1,1)) -(line3(2:end,1)-line3(1:end-1,1));

elast =((line1(2:end,end)-line3(2:end,end))-(line1(1:end-1,end)-line3(1:end-1,end)))./(line1(2:end,5)-line1(1:end-1,5));
test1 = ((line1(2:end,end)-line3(2:end,end))-(line1(1:end-1,end)-line3(1:end-1,end)));
test2 = (line1(2:end,5)-line1(1:end-1,5))*100;
test = [test1, test2]
figure
plot(test,'DisplayName','test')