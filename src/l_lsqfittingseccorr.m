%% The same thing as l_lsqfitting but with the second correction
%CMM 2020



%% For only phi vs H
a0 = [1 7];

H1=datastack.H;
plotthison=datastack.H;%or H1
%set(0,'defaultAxesFontSize',10)

%
wherehighGBD = datastack.GBDseccorr>7;
H1=datastack.H.*wherehighGBD;
H1(H1==0)=NaN;
%}
%turn into column of data and remove nan
xdata = datastack.Phireflseccorr(:)*180/pi();
H2fit=H1(:);
xdata=xdata(~isnan(H2fit));
xlocleft=X(~isnan(H2fit)); %keeping the x and y for stuff
ylocleft=Y(~isnan(H2fit));
H2fit=H2fit(~isnan(H2fit));

[a, resnorm,residual,~,~,~,J] = lsqcurvefit(@myfun,a0,xdata,H2fit); %least square
hrsq=1-sum(residual.^2)/sum((H2fit-mean(H2fit)).^2); %1-sum squared regression/sum of squares

[ypred,delta] = nlpredci(@myfun,xdata,a,residual,'Jacobian',J,'PredOpt','observation'); %do a prediction
lower = ypred-delta;
upper = ypred+delta;
dataupperlower=[xdata,upper,lower];

figure;
data2=linspace(10,90);
F=myfun(a,data2);
%F=a(1)*cos((pi()/90)*data2(:))+a(2);
plot(data2,F,'r','LineWidth',2)
hold on
dataupperlower=sortrows(dataupperlower);
plot(dataupperlower(:,1),[dataupperlower(:,2),dataupperlower(:,3)],'r:', 'linewidth',2)
scatter(xdata(:),H2fit(:),'bx')
txt = (['H = ',num2str(a(1),'%4.2f'), '*cos(\theta)+', num2str(a(2),'%4.2f')]); %
text(55, 4.2, txt);%, 'FontSize', 10, 'Color', 'k');
txt2 = (['R^{2} = ',num2str(hrsq,'%4.2f')]); %
text(55, 4.0, txt2);%, 'FontSize', 10, 'Color', 'k')
legend({'Cosine Fit','95% confidence interval'},'Location','northeast')
hold off 
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5.5*nanstd(H1(:))])
xlim([0 90])
title('Declination angle against measured hardness with fit')
figname=['lowGBD Phi V H Figure FlippedFITCOS SECOND CORR_ALL' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end
%Basics only: comment out:
%{

%% Looking at outliers: where are they?


%datacolumn=[xdata,H2fit,xlocleft,ylocleft];
datacolumn=[datastack.Phireflseccorr(:)*180/pi(),datastack.H(:),datastack.X(:),datastack.Y(:)];

%now data is sorted in an easy way to compare to the function plus or minus
%the mean delta
%first find points that are too high
locwheretoohigh=[];
for i=1:size(datacolumn,1)
    if datacolumn(i,2)>(mean(delta)+myfun(a,datacolumn(i,1)))
        locwheretoohigh(size(locwheretoohigh,1)+1,1:2)= datacolumn(i,3:4);
    end %if it's bigger, pull it out and save the location
end
figoutHhigh=figure;
figoutHhigh.Name='outHhigh';
hplot=contourf(X,Y,plotthison,45,'LineColor','None');
hold on
title('Nanoindentation map highlighting points beyond 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
scatter(locwheretoohigh(:,1),locwheretoohigh(:,2),'rx');
hold off
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['outHhigh2c' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end

%now find points too low
locwheretoolow=[];
for i=1:size(datacolumn,1)
    if datacolumn(i,2)<(-mean(delta)+myfun(a,datacolumn(i,1)))
        locwheretoolow(size(locwheretoolow,1)+1,1:2)= datacolumn(i,3:4);
    end
end
figoutHlow=figure;
figoutHlow.Name='outHlow';
hplot=contourf(X,Y,plotthison,45,'LineColor','None');
hold on
title('Nanoindentation map highlighting points below 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
scatter(locwheretoolow(:,1),locwheretoolow(:,2),'kx');
hold off
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['outHlow2c' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end

%where are these points w.r.t. the phi? 
figoutPHIlow=figure;
figoutPHIlow.Name='outPHIlow';
hplot=contourf(X,Y,datastack.Phiseccorr*180/pi,45,'LineColor','None');
hold on
title('EBSD map highlighting points above and below 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
ha=scatter(locwheretoohigh(:,1),locwheretoohigh(:,2),'rx');
hb=scatter(locwheretoolow(:,1),locwheretoolow(:,2),'kx');
legend(([ha hb]),'Points above 95% CI','Points below 95% CI')
hold off
c=colorbar;
c.Label.String = 'Phi angle /deg';
figname=['outPHIlow2c' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end
%% variation of grain boundary distance filters effect on R squared of fit

hrsqvariablegbd=zeros(20,5);
round(max(max(datastack.GBDseccorr)))
for i=1:20
    wherehighGBDVAR = datastack.GBDseccorr>(i-1)*round(max(max(datastack.GBDseccorr)))/20;
    H1var=datastack.H.*wherehighGBDVAR;
    H1var(H1var==0)=NaN;
    xdatavar = datastack.Phireflseccorr(:)*180/pi;
    H2fitvar=H1var(:);
    xdatavar=xdatavar(~isnan(H2fitvar));
    H2fitvar=H2fitvar(~isnan(H2fitvar));
    [avar,~,residualvar] = lsqcurvefit(@myfun,a0,xdatavar,H2fitvar); %least square
    hrsqvariablegbd(i,2)=1-sum(residualvar.^2)/sum((H2fitvar-mean(H2fitvar)).^2); %1-sum squared regression/sum of squares
    hrsqvariablegbd(i,1)=(i-1)*round(max(max(datastack.GBDseccorr)))/20;
    hrsqvariablegbd(i,3:4)=avar;    
end

figphiVHlgbd=figure;
figphiVHlgbd.Name='varGBD vs Rsquared';
scatter(hrsqvariablegbd(:,1),hrsqvariablegbd(:,2))
ylim([0 1])
xlabel('Exclusion zone size /\mum')
ylabel('R squared of fit')
title({'Effect of grain boundary distance filter';'on R squared of fit'})
figname=['varGBD V Rsquared Figure 2c' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end
%}
close all
%% functions

function F = myfun(a,data)
    x = data(:);
    F=a(1)*cos((pi()/90)*x)+a(2);
end