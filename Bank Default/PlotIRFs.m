clear all
close all
clc
addpath('c:\dynare\4.5.7\matlab');
addpath('C:\Users\Farkas\Dropbox\PhD\Research\QE effects\Full calibration\export_fig-master');
dynare RP_omega



%% Technology shock
figure
subplot(3,3,1)
plot([1:40],oo_.irfs.Y_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('$Y_t$, Output','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,2)
plot([1:40],oo_.irfs.I_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('$I_t$, Investment','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,3)
plot([1:40],oo_.irfs.q_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('$q_t$, Price of Capital','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,4)
plot([1:40],oo_.irfs.bigN_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('$N_t$, Entrepreneurial Net Worth','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,5)
plot([1:40],oo_.irfs.mu_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('')
title('$\mu_t$ Monitoring Intensity','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,6)
plot([1:40],oo_.irfs.BVaR_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('')
title('$\Gamma_t$, Bank Default Probability','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')



subplot(3,3,7)
plot([1:40],oo_.irfs.TL_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('$I_t-N_t$, Total Lending','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,8)
plot([1:40],oo_.irfs.bigA_z_shk/sigmaz,'color',[0 0 0],'LineWidth',1);
title('')
title('$A_t$, Bank Capital','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,9)
plot([1:40],1./(oo_.irfs.gammag_z_shk/sigmaz+gammagtemp)*100,'color',[0 0 0],'LineWidth',1);
title('$\gamma_t$, Regulatory Capital Ratio','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, Level')
xlabel('Quarters')

sgtitle('Impulse Responses to a Negative Technology Shock in Normal Times');
sgt.FontSize = 14; sgt.FontWeight = 'Bold';
set(gcf,'color','w');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

export_fig('TechnologyShock.pdf')


%%
figure
subplot(3,3,1)
plot([1:40],oo_.irfs.Y_bk_shk,'color',[0 0 0],'LineWidth',1);
title('$Y_t$, Output','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,2)
plot([1:40],oo_.irfs.I_bk_shk,'color',[0 0 0],'LineWidth',1);
title('$I_t$, Investment','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,3)
plot([1:40],oo_.irfs.q_bk_shk,'color',[0 0 0],'LineWidth',1);
title('$q_t$, Price of Capital','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,4)
plot([1:40],oo_.irfs.bigA_bk_shk,'color',[0 0 0],'LineWidth',1);
title('$A_t$, Bank Capital','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,5)
plot([1:40],oo_.irfs.mu_bk_shk,'color',[0 0 0],'LineWidth',1);
title('')
title('$\mu_t$, Monitoring Intensity','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,6)
plot([1:40],oo_.irfs.BVaR_bk_shk,'color',[0 0 0],'LineWidth',1);
title('')
title('$\Gamma_t$, Bank Default Probability','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')



subplot(3,3,7)
plot([1:40],oo_.irfs.TL_bk_shk,'color',[0 0 0],'LineWidth',1);
title('$I_t-N_t$, Total Lending','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,3,8)
plot([1:40],oo_.irfs.Rd_bk_shk,'color',[0 0 0],'LineWidth',1);
title('')
title('$R_d$, Policy Rate','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

subplot(3,3,9)
plot([1:40],1./(oo_.irfs.gammag_bk_shk+gammagtemp)*100,'color',[0 0 0],'LineWidth',1);
title('$\gamma_t$, Regulatory Capital Ratio','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, Level')
xlabel('Quarters')


sgtitle('Impulse Responses to a Negative Bank Capital Shock in Normal Times');
sgt.FontSize = 14; sgt.FontWeight = 'Bold';
set(gcf,'color','w');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

export_fig('BankCapitalShock.pdf')

%%
figure

subplot(2,2,1)

plot([1:40],-oo_.irfs.Rd_mp_shk,'color',[0 0 0],'LineWidth',1);
title('$a_t$, Bank Capital','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(2,2,2)

plot([1:40],-oo_.irfs.BVaR_mp_shk,'color',[0 0 0],'LineWidth',1);
title('$\Gamma_t$, Bank Default Probability','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')



subplot(2,2,3)
plot([1:40],-oo_.irfs.mu_mp_shk,'color',[0 0 0],'LineWidth',1);
title('$\mu_t$, Monitoring Intesity','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(2,2,4)
plot([1:40],-oo_.irfs.q_mp_shk,'color',[0 0 0],'LineWidth',1);
title('')
title('$q_t$, Real Price of Capital','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')

sgtitle('Impulse Responses to a Expansionary Monetary Policy Shock in Normal Times');
sgt.FontSize = 14; sgt.FontWeight = 'Bold';
set(gcf,'color','w');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

export_fig('MonPolShock.pdf')
%%
%%
figure

subplot(3,2,1)

plot([1:40],oo_.irfs.gammag_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$\gamma_t$, Capital Requirement','interpreter','latex')
% ylim([-0.; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,2,2)

plot([1:40],oo_.irfs.BVaR_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$\Gamma_t$, Bank Default Probability','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')



subplot(3,2,3)
plot([1:40],oo_.irfs.mu_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$\mu_t$, Monitoring Intesity','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,2,4)
plot([1:40],oo_.irfs.q_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$q_t$, Real Price of Capital','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,2,5)
plot([1:40],oo_.irfs.Rd_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$R_d$, Policy Rate','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')


subplot(3,2,6)
plot([1:40],oo_.irfs.TL_macropru_shk,'color',[0 0 0],'LineWidth',1);
title('$I_t-N_t$, Total Lending','interpreter','latex')
% ylim([0; 0.25])
ylabel('%, deviation from SS')
xlabel('Quarters')



sgtitle('Impulse Responses to a Macro Prudential Policy Easing Shock in Normal Times');
sgt.FontSize = 14; sgt.FontWeight = 'Bold';
set(gcf,'color','w');
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');

export_fig('MacroPruShock.pdf')