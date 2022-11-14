function paperpplot1(titlelist,legendlist,figlabel,ylabels,zdata1,zdata2,zdata3)



figure
set(gcf,'color','w');

% titlelist = char(strrep(cellstr(titlelist),'_','.'));

ndsets=3;       % default, changed below as applicable
if nargin==5
    zdata2=nan*zdata1;
    zdata3=nan*zdata1;
    ndsets =1;
elseif nargin == 6
    zdata3 =nan*zdata1;
    ndsets=2;
elseif ((nargin>8) | (nargin <=4))
    error ('makechart takes 5 to 6 arguments')
end

nobs = size(zdata1,1);
xvalues = (1:nobs)';

nvars = size(titlelist,1);
if nvars==1 
    nrows=1;
    ncols = 1;
elseif nvars==2
    nrows =2;
    ncols = 1;
elseif (nvars == 3 | nvars ==4)
    nrows = 2;
    ncols =2;
elseif (nvars==5 |nvars ==6)
    nrows = 3;
    ncols = 2;
elseif (nvars==7 | nvars==8)
    nrows = 4;
    ncols = 2;
elseif (nvars==9 | nvars==10)
    nrows = 5;
    ncols = 2;
else 
    error('too many variables (makechart)')
end

for i = 1:nvars
    subplot(nrows,ncols,i)
    h1=plot(xvalues,zdata1(:,i),'r-','linewidth',2); hold on
    h1=plot(xvalues,zdata2(:,i),'k--','linewidth',2); hold on
    h2=plot(xvalues,zdata3(:,i),'k-','LineWidth',2);
    [x0 x1 y10 y11] = pickaxes(xvalues,zdata1(:,i));
    [x0 x1 y20 y21] = pickaxes(xvalues,zdata2(:,i));
    [x0 x1 y30 y31] = pickaxes(xvalues,zdata3(:,i));
    y0 = min([y10,y20,y30]);
    y1 = max([y11,y21,y31]);
    if y0==y1
        y1=y0+1;
    end
    
    axis([x0 x1 y0 y1])
    set(h1);
 
    if i==nvars | i==nvars-1
        xlabel('Time');
    end
%     set(gca,'XTick',xtick)
%     set(gca,'XTickLabel',xticklabel)
    
    title([num2str(i),'. ',titlelist(i,:)],'interpreter','latex');
    ylabel(ylabels(i,:));
  
end
      hL = legend(legendlist);
% Programatically move the Legend
newPosition = [0.5 0.05 0 0];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits, 'Orientation','horizontal');
% sets printing preferences
%printpref
