%%Function to run some least squares fitting to data on low oxygen hardness
%%values (H1)
%CMM 2019


%% For oxygen data

a0 = [1 1];


%remove saturated titanium
%XPCcontourf(xdata)
%{
wheremat = datastack.EPMAO > 10;
Hmatunsat=Hmat.*~wheremat;
Hmatunsat(Hmatunsat==0)=NaN;
datastack.Hmatunsat=Hmatunsat;
%}
datastack.Hmatunsat=datastack.Hmat;
%turn into column of data and remove nan
xdata = double(datastack.EPMAOS); %set to xdata 
%normalise chemistry based on average O level in first 5 columns
xshift=mean(mean(xdata(1:2,:)));
xdata=xdata-xshift;
xdata2=zeros(size(xdata,1)*size(xdata,2),2);
xdata2(:,1) = xdata(:);
xdata2(:,2) = datastack.Phirefl(:)*180/pi;


H2fit = datastack.Hmatunsat(:);

xdata2=xdata2(~isnan(H2fit),:);
H2fit=H2fit(~isnan(H2fit));

%quick remove of x<0 points
xdata2(xdata2(:,1)<0)=NaN;
H2fit=H2fit(~isnan(xdata2(:,1)));
xdata2=xdata2(~isnan(xdata2(:,1)),:);

%remove points between x=8 and 10
xdata2(xdata2(:,1)>7.5 & xdata2(:,1)<11)=NaN;
H2fit=H2fit(~isnan(xdata2(:,1)));
xdata2=xdata2(~isnan(xdata2(:,1)),:);



[a, resnorm,residual,~,~,~,J] = lsqcurvefit(@oxygenfitwORI,a0,xdata2,H2fit); %least square
hrsq=1-sum(residual.^2)/sum((H2fit-mean(H2fit)).^2); %1-sum squared regression/sum of squares

J=full(J);

[ypred,delta] = nlpredci(@oxygenfitwORI,xdata2,a,residual,'Jacobian',J,'PredOpt','observation'); %do a prediction
lower = ypred-delta;
upper = ypred+delta;
dataupperlower=[xdata2,upper,lower];


%set up the model surface
[data2x,data2y] = meshgrid(min(xdata2(:,1)):range(xdata2(:,1))/100:max(xdata2(:,1)),min(xdata2(:,2)):range(xdata2(:,2))/100:max(xdata2(:,2)));
%[data2x,data2y] = meshgrid(min(xdata2(:,1)):range(xdata2(:,1))/100:max(xdata2(:,1)),...
%    0:90/100:90);
data2=[data2x(:) data2y(:)];
F=oxygenfitwORI(a,data2);
F=gridify_vector(F',size(data2x,1),size(data2x,1));
F=F'; %flip because things

%% reevaluate x and y data for plot

datastack.Hmatunsat=datastack.Hmat;
%turn into column of data and remove nan
xdata = double(datastack.EPMAOS); %set to xdata 
%normalise chemistry based on average O level in first 5 columns
xshift=mean(mean(xdata(1:2,:)));
xdata=xdata-xshift;
xdata2=zeros(size(xdata,1)*size(xdata,2),2);
xdata2(:,1) = xdata(:);
xdata2(:,2) = datastack.Phirefl(:)*180/pi;


H2fit = datastack.Hmatunsat(:);

xdata2=xdata2(~isnan(H2fit),:);
H2fit=H2fit(~isnan(H2fit));

%quick remove of x<0 points
xdata2(xdata2(:,1)<0)=NaN;
H2fit=H2fit(~isnan(xdata2(:,1)));
xdata2=xdata2(~isnan(xdata2(:,1)),:);
%% plot

figure;
surf(data2x,data2y,F,data2y,'FaceAlpha',0.7,'EdgeColor','None')
c=colorbar;
c.Label.String = 'Declination angle /^{o}';
hold on
scatter3(xdata2(:,1),xdata2(:,2),H2fit(:),20,xdata2(:,2),'filled','MarkerEdgeColor','k')
txt = (['H = ',num2str(a(1),'%4.2f'), '* O^{' ,num2str(a(2),'%4.2f') ,...
    '}+0.67*cos(\theta)+3.2']);
text(11, 50,5, txt);%, 'FontSize', 10, 'Color', 'k');
txt2 = (['R^{2} = ',num2str(hrsq,'%4.2f')]); %
text(11, 50,3, txt2);%, 'FontSize', 10, 'Color', 'k');
hold off 
xlabel('Oxygen Concentration / arb units')
ylabel('Declination Angle /^{o}')
zlabel('Hardness /GPa')
xlim([0 30])
ylim([50 90])
zlim([0 35])
title('Oxygen and declination angle against measured hardness with fit')
figname=['OX V H Figure with Ori' ebsdname(1:(max(size(ebsdname)-4)))];
%print(fullfile(resultsdir, figname),'-dpng',resolution)

%saveas(gcf,fullfile(resultsdir, figname),'png')
%saveas(gcf,fullfile(resultsdir, figname),'fig')
 nnz(~isnan(H2fit(:))) %count how many are in the fig


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
%close all
%XPCcontourf(datastack.H)
%% functions

function F = myfun(a,data)
    x = data(:);
%    F=a(1)*x.^3+a(2)*x.^2+a(3)*x+a(4);
    F=a(1)*((x).^a(2));
end