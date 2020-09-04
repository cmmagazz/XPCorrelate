%%Function to run some least squares fitting to data on low oxygen hardness
%%values (H1)
%CMM 2019


%% For oxygen data

a0 = [1 1 0];


%remove saturated titanium
%{
wheremat = datastack.EPMAO > 10;
Hmatunsat=Hmat.*~wheremat;
Hmatunsat(Hmatunsat==0)=NaN;
datastack.Hmatunsat=Hmatunsat;
%}
datastack.Hmatunsat=datastack.Hmat;
%turn into column of data and remove nan
xdata = double(datastack.EPMAO(:));
H2fit = datastack.Hmatunsat(:);
xdata=xdata(~isnan(H2fit));
xlocleft=datastack.X(~isnan(H2fit)); %keeping the x and y for stuff
ylocleft=datastack.Y(~isnan(H2fit));
H2fit=H2fit(~isnan(H2fit));

[a, resnorm,residual,~,~,~,J] = lsqcurvefit(@myfun,a0,xdata,H2fit); %least square
hrsq=1-sum(residual.^2)/sum((H2fit-mean(H2fit)).^2); %1-sum squared regression/sum of squares

J=full(J);

[ypred,delta] = nlpredci(@myfun,xdata,a,residual,'Jacobian',J,'PredOpt','observation'); %do a prediction
lower = ypred-delta;
upper = ypred+delta;
dataupperlower=[xdata,upper,lower];

figure;
data2=linspace(min(xdata),max(xdata));
F=myfun(a,data2);
plot(data2,F,'r','LineWidth',2)
hold on
dataupperlower=sortrows(dataupperlower);
plot(dataupperlower(:,1),[dataupperlower(:,2),dataupperlower(:,3)],'r:', 'linewidth',2)
scatter(xdata(:),H2fit(:),'bo')
txt = (['H = ',num2str(a(1),'%4.2f'), '* \surd{O} + ' ,num2str(a(2),'%4.2f') ,'*O + ', num2str(a(3),'%4.2f')]);
text(11, 2, txt);%, 'FontSize', 10, 'Color', 'k');
legend({'Fit','95% confidence interval'},'Location','northeast')
hold off 
xlabel('Oxygen Concentration / arb units')
ylabel('Hardness /GPa')
ylim([0 35])
xlim([0 30])
title('Oxygen against measured hardness with fit')
figname=['OX V H Figure FlippedFITCOS' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')


%% Looking at outliers: where are they?
%{

datacolumn=[xdata,H2fit,xlocleft,ylocleft];
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
hplot=contourf(datastack.X,datastack.Y,datastack.Hmat,45,'LineColor','None');
hold on
title('Nanoindentation map highlighting points beyond 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
scatter(locwheretoohigh(:,1),locwheretoohigh(:,2),'rx');
hold off
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['outHhigh' ebsdname(1:(max(size(ebsdname)-4)))];
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
hplot=contourf(datastack.X,datastack.Y,datastack.Hmat,45,'LineColor','None');
hold on
title('Nanoindentation map highlighting points below 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
scatter(locwheretoolow(:,1),locwheretoolow(:,2),'kx');
hold off
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['outHlow' ebsdname(1:(max(size(ebsdname)-4)))];
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
figoutgIDlow=figure;
figoutgIDlow.Name='outgIDlow';
hplot=contourf(datastack.X,datastack.Y,datastack.gID,max(max(datastack.gID)),'LineColor','None');
hold on
title('Nanoindentation map highlighting points below 95% CI')
xlabel('\mum')
ylabel('\mum')
axis image
scatter(locwheretoolow(:,1),locwheretoolow(:,2),'kx');
hold off
c=colorbar;
c.Label.String = 'Grain ID';
figname=['outgIDlow' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end
%}
close all
XPCcontourf(datastack.H)
%% functions

function F = myfun(a,data)
    x = data(:);
    F=a(1)*x.^0.5+a(2)*x+a(3);
end