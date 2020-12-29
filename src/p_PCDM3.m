%% PCDM3
%Script with even more PCDM, in response to requests from reviewers (Oct
%2020) 
%
% CMM 2020
%
%% 2D Histogram plot of H vs Phi
%
hist3([datastack.Phireflseccorr(:).*180/pi(), datastack.H(:)],'CdataMode','auto','nbins',[30 30])
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([2 4.8])
xlim([15 90])
c=colorbar;
c.Label.String='Number of indents';
view(2)
title('2D Histogram of declination angle against measured hardness')
figname=['2D histogram phi vs h' filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
%% Extracting a specific grain
%
grainno=90;
setofinterest =  datastack.gIDseccorr == grainno;

Hofinterest=datastack.H.*setofinterest;
Hofinterest(Hofinterest==0)=NaN;

XPCcontourf(datastack.gIDseccorr,'title',"gID")

XPCcontourf(Hofinterest,'title',"Grain 69")

figure;
hist=histogram(Hofinterest(:));
title(['Histogram of Hardness Measurements in Grain ' num2str(grainno)])
xlabel('Hardness /GPa')
xlim([min(Hofinterest(:)) max(Hofinterest(:))]) 
ylabel('Number of Indents')
txt = {['Average Hardness: ' num2str(nanmean(Hofinterest(:)), '%.3g') ' GPa'],...
    ['Standard Deviation: ' num2str(nanstd(Hofinterest(:)), '%.3g') ' GPa']};
text(0.05*max(Hofinterest(:)),max(hist.Values(:))*0.9,txt)
figname=['Hardness Histogram in grain ' num2str(grainno) filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)

%% Comparisons of above
XPCcontourf(datastack.gIDseccorr)

setofinterest =  datastack.gIDseccorr == 34;
Hofinterest34=datastack.H.*setofinterest;
Hofinterest34(Hofinterest34==0)=NaN;

phiofinterest34=datastack.Phireflseccorr.*setofinterest;
phiofinterest34(phiofinterest34==0)=NaN;
nanmean(phiofinterest34(:))*180/pi

setofinterest =  datastack.gIDseccorr == 63;
Hofinterest63=datastack.H.*setofinterest;
Hofinterest63(Hofinterest63==0)=NaN;

phiofinterest63=datastack.Phireflseccorr.*setofinterest;
phiofinterest63(phiofinterest63==0)=NaN;
nanmean(phiofinterest63(:))*180/pi

figure;
hist1=histogram(Hofinterest34(:));
hold on
hist2=histogram(Hofinterest63(:));
hist1.BinWidth = 0.075;
hist2.BinWidth = 0.075;
title('Histogram of Hardness Measurements in Grains 34 and 63')
xlabel('Hardness /GPa')
txt = {['34: Average Hardness: ' num2str(nanmean(Hofinterest34(:)), '%.3g') ' GPa'],...
    ['Standard Deviation: ' num2str(nanstd(Hofinterest34(:)), '%.3g') ' GPa']};
text(0.4*max(Hofinterest34(:)),max(hist1.Values(:))*0.9,txt)
txt2 = {['63: Average Hardness: ' num2str(nanmean(Hofinterest63(:)), '%.3g') ' GPa'],...
    ['Standard Deviation: ' num2str(nanstd(Hofinterest63(:)), '%.3g') ' GPa']};
text(0.4*max(Hofinterest63(:)),max(hist2.Values(:))*0.9,txt2)
%xlim([hist1.BinLimits(1) hist1.BinLimits(2)])
ylabel('Number of Indents')
legend({'Grain 34', 'Grain 63'});
figname=['Hardness Histogram in of grains 34 and 63' filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)

%% Modulus
%}
setofinterest =  datastack.gIDseccorr == 34;
Mofinterest34=datastack.M.*setofinterest;
Mofinterest34(Mofinterest34==0)=NaN;

setofinterest =  datastack.gIDseccorr == 63;
Mofinterest63=datastack.M.*setofinterest;
Mofinterest63(Mofinterest63==0)=NaN;

figure;
hist1=histogram(Mofinterest34(:));
hold on
hist2=histogram(Mofinterest63(:));
hist1.BinWidth = 4;
hist2.BinWidth = 4;
title('Histogram of Modulus Measurements in Grains 34 and 63')
xlabel('Modulus /GPa')
txt = {['34: Average Modulus: ' num2str(nanmean(Mofinterest34(:)), '%.3g') ' GPa'],...
    ['Standard Deviation: ' num2str(nanstd(Mofinterest34(:)), '%.3g') ' GPa']};
text(0.5*max(Mofinterest34(:)),max(hist1.Values(:))*0.9,txt)
txt2 = {['63: Average Modulus: ' num2str(nanmean(Mofinterest63(:)), '%.3g') ' GPa'],...
    ['Standard Deviation: ' num2str(nanstd(Mofinterest63(:)), '%.3g') ' GPa']};
text(0.5*max(Mofinterest63(:)),max(hist2.Values(:))*0.9,txt2)
%xlim([hist1.BinLimits(1) hist1.BinLimits(2)])
ylabel('Number of Indents')
legend({'Grain 34', 'Grain 63'});
figname=['Modulus histogram of grains 34 and 63' filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)


%% Isolating grainno 69 and 90 and plotting those on the second 
% corrected scatter plot with different colours
figure;
scatter(datastack.Phireflseccorr(:)*180/pi(),datastack.H(:),'x')
hold on
scatter(datastack.Phireflseccorr(:)*180/pi(),Hofinterest69(:),'x')
scatter(datastack.Phireflseccorr(:)*180/pi(),Hofinterest90(:),'x')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(datastack.H(:))-5*nanstd(datastack.H(:)) nanmean(datastack.H(:))+5*nanstd(datastack.H(:))])
xlim([0 90])
title('Declination angle against measured hardness')
legend({'Full map','Grain 69', 'Grain 90'})
figname=['SECCORR phi vs H grain69 and 90' ebsdname(1:(max(size(ebsdname)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
%% Plotting histograms of similar orientations vs 69 and 90 respectively
%first 69 
setofinterest =  datastack.gIDseccorr == 69;
Hofinterest69=datastack.H.*setofinterest;
Hofinterest69(Hofinterest69==0)=NaN;

setofinterest =  datastack.gIDseccorr ~= 69 & datastack.Phireflseccorr> 26*pi/180 & datastack.Phireflseccorr< 36*pi/180;
Hofinterest69neigh=datastack.H.*setofinterest;
Hofinterest69neigh(Hofinterest69neigh==0)=NaN;

XPCcontourf(Hofinterest69neigh,'title',"Hardness of Grains within 5^{o} of Grain 69",...
    'saveq',1,'cunits',"Hardness /GPa")

figure;
hist1=histogram(Hofinterest69(:));
hold on
hist2=histogram(Hofinterest69neigh(:));
hist1.BinWidth = 0.1;
hist2.BinWidth = 0.1;
title({'Histogram of Hardness Measurements in Grains 69',' and grains within 5^{o} of Grain 69'})
xlabel('Hardness /GPa')
%xlim([hist1.BinLimits(1) hist1.BinLimits(2)])
ylabel('Number of Indents')
legend({'Grain 69', 'Grains within 5^{o} of Grain 69'});
figname=['Hardness Histogram in of grains 69 and 69neigh' filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)

%% Now 90
setofinterest =  datastack.gIDseccorr == 90;
Hofinterest90=datastack.H.*setofinterest;
Hofinterest90(Hofinterest90==0)=NaN;

setofinterest =  datastack.gIDseccorr ~= 90 & datastack.Phireflseccorr> 51*pi/180 & datastack.Phireflseccorr< 61*pi/180;
Hofinterest90neigh=datastack.H.*setofinterest;
Hofinterest90neigh(Hofinterest90neigh==0)=NaN;

XPCcontourf(Hofinterest90neigh,'title',"Hardness of Grains within 5^{o} of Grain 69",...
    'saveq',1,'cunits',"Hardness /GPa")

figure;
hist1=histogram(Hofinterest90(:));
hold on
hist2=histogram(Hofinterest90neigh(:));
hist1.BinWidth = 0.1;
hist2.BinWidth = 0.1;
title({'Histogram of Hardness Measurements in Grains 90',' and grains within 5^{o} of Grain 90'})
xlabel('Hardness /GPa')
%xlim([hist1.BinLimits(1) hist1.BinLimits(2)])
ylabel('Number of Indents')
legend({'Grain 90', 'Grains within 5^{o} of Grain 90'});
figname=['Hardness Histogram in of grains 90 and neigh90' filename(1:(max(size(filename)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
%}

%% Looking at the Fused silica map
uiopen('Z:\CM\18_OctEXPRESS\191203_ExpressFUSI\indentdataonly.xlsx',1)
%indentation data is all in a big column based on order that the indent was
%performed. Modulus is in column 4, hardness in 6. 
windowwidth=200;
indentno=linspace(1,size(indentdataonly,1),size(indentdataonly,1))';
smoothM=smooth(indentdataonly(:,4),windowwidth);
smoothH=smooth(indentdataonly(:,6),windowwidth);
Mstd=movstd(indentdataonly(:,4),windowwidth);
Hstd=movstd(indentdataonly(:,6),windowwidth);

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)
hplot=shadedErrorBar(indentno(:),smoothH(:),Hstd(:),'lineProps',{'b-','markerfacecolor','b'})
xlabel('Indent Number')
ylabel('Hardness (GPa)')
hold off


subplot(1,2,2)
hplot=shadedErrorBar(indentno(:),smoothM(:),Mstd(:),'lineProps',{'b-','markerfacecolor','b'})
xlabel('Indent Number')
ylabel('Modulus (GPa)')
hold off

sgtitle({'Scatter plot of nanoindentation map in fused silica,','using a smoothing window of 200 indents'})

figname='fusimap_trend';
resolution = ['-r' num2str(1500)];
print(fullfile(resultsdir, figname),'-dpng',resolution)

%% Line profiles across EPMA and Hardness maps
wheremat = datastack.H < 12 & datastack.X > 160;
Hmat2=H.*~wheremat;
wheremat = datastack.H < 5 & datastack.X > 150;
Hmat2=Hmat2.*~wheremat;
EPMAOSshift=datastack.EPMAOSshift;
X=datastack.X;

Hmat2(Hmat2==0)=NaN;
EPMAOSshift(isnan(Hmat2))=NaN;
X(isnan(Hmat2))=NaN;
Xave=nanmean(X,2);

EPMAaverage=nanmean(datastack.EPMAOSshift,2);
EPMAstd=nanstd(datastack.EPMAOSshift,0,2);

Haverage=nanmean(Hmat2,2);
Hstd=nanstd(Hmat2,0,2);

figure()
yyaxis left
shadedErrorBar(Xave(:),EPMAaverage(:),EPMAstd(:),'lineProps',{'b-','markerfacecolor','b'})
ylim([0,18])
ylabel('Oxygen Concentration /arb units')
hold on
yyaxis right
ylabel('Hardness /GPa')
shadedErrorBar(Xave(:),Haverage(:),Hstd(:),'lineProps',{'r-','markerfacecolor','r'})
xlabel('Distance from edge /\mum')
title({'Integrated line profiles along EPMA' 'and Nanoindentation hardness map'})

figname='lineintprofile_epmavsO';
resolution = ['-r' num2str(1500)];
print(fullfile(resultsdir, figname),'-dpng',resolution)