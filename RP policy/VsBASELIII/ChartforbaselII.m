clear all
close all
clc
addpath('c:\dynare\4.5.7\matlab');
addpath('C:\Users\Farkas\Dropbox\PhD\Research\QE effects\Full calibration\export_fig-master');


B3u = 0.08*ones(1,11)+ [linspace(0,0.025,9) 0.025 0.025];
B3d =0.08*ones(1,11) - [linspace(0,0.025,9) 0.025 0.025];
B3 = [ fliplr(B3d) 0.08*ones(1,3)  B3u]
xaxis = linspace(-0.12,0.12,25);
myrule=1*ones(1,25)./(12.5*ones(1,25)-32.46*xaxis)

plot(xaxis,B3,'--k','Linewidth',1)
axis tight
hold on 
plot(xaxis,myrule,'-k','Linewidth',2)

     a=[cellstr(num2str(get(gca,'xtick')'*100))]; 
% Create a vector of '%' signs
     pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
     new_xticks = [char(a),pct];
% 'Reflect the changes on the plot
     set(gca,'xticklabel',new_xticks)
     % Convert y-axis values to percentage values by multiplication
     a=[cellstr(num2str(get(gca,'ytick')'*100))]; 
% Create a vector of '%' signs
     pct = char(ones(size(a,1),1)*'%'); 
% Append the '%' signs after the percentage values
     new_yticks = [char(a),pct];
% 'Reflect the changes on the plot
     set(gca,'yticklabel',new_yticks)
     xlabel('Credit-to-GDP gap','FontSize',16);
     ylabel('Regulatory Capital Ratio','FontSize',16);
     legend('Basel III', 'Model Macroprudential Rule', 'Location','southeast','FontSize',16);
     title('Macroprudential Policy Regulation','FontSize',16)

     
                mkdir('Figures');
                cd Figures
                
                frame_h = get(handle(gcf),'JavaFrame');
                set(frame_h,'Maximized',1);
                set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
                % Get rid of tool bar and pulldown menus that are along top of figure.
                set(gcf, 'Toolbar', 'none', 'Menu', 'none');

                saveas(gcf,['FigureBASELIII.fig']);
                print('FigureBASELIII','-dpng');
                matlab2tikz('BaselIII.tex');