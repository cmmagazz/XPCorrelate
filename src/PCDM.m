%A collection of bits of code to play around with data once you've
%corrected it. I.e. Post-Corrected Data Manipulator.
%CMM 2019
%looking at low oxygen (therefore bulk) relationships of phi vs H

%first fix the phi to be between 0 and 90 FOR HCP?
datastack.Phi(datastack.Phi>(pi/2))=pi-datastack.Phi(datastack.Phi>(pi/2));

%% Remove points close to a grain boundary from phi vs H plot 
wherehighGBD = datastack.GBD>7;
H1=datastack.H.*wherehighGBD;
H1(H1==0)=NaN;
wherehighGBD = datastack.GBD<2;
H2=datastack.H.*wherehighGBD;
H2(H2==0)=NaN;

%number of points
sum(wherehighGBD(:))

figphiVHlgbd=figure;
figphiVHlgbd.Name='Phi vs H';
scatter(datastack.Phirefl(:)*180/pi(),H1(:) ,'x')
xlabel('Declination angle  /^{o}')
ylabel('Hardness /GPa')
ylim([1 5])
xlim([0 90])
title({'Declination angle against measured nanoindentation hardness,';' excluding points near grain boundaries'})
figname=['LowGBD Phi V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
Hlgbd=H1;
XPCcontourf(Hlgbd,'title',"Nanoindentation hardness - excluding points near GB",...
    'cunits',"Hardness /GPa",'saveq',1)


figphiVHlgbd=figure;
figphiVHlgbd.Name='Phi vs H';
scatter(datastack.Phirefl(:)*180/pi,H1(:))
hold on
scatter(datastack.Phirefl(:)*180/pi,H2(:),'rx')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
legend({'Points closer than 7um (inverse of fig11)','Points closer than 3um'})
ylim([1.5 5])
title({'Declination angle against measured hardness';'displaying only points close to grain boundaries'})
figname=['lowGBD Phi V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')
XPCcontourf(H1)
XPCcontourf(H2)

%{
whereLowGBD = GBD<7;
Hgb=H.*whereLowGBD;
Hgb(Hgb==0)=NaN;

figure;
contourf(X,Y,Hmat,45,'LineColor','None');
title('H Near Grain Boundaries')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['HNGB' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
%}
figH=figure;
figH.Name='HABG';
contourf(X,Y,H1,45,'LineColor','None');
title('Nanoindentation map excluding points near grain boundaries')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['HABG' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

whereLowGBD = GBD>7;
M1=M.*whereLowGBD;
M1(M1==0)=NaN;
figphiVHlgbd=figure;
figphiVHlgbd.Name='Phi vs M';
scatter(FPhiGrtest(:)*180/pi,M1(:))
xlabel('Phi /degrees')
ylabel('Modulus /GPa')
%ylim([0 6])
title('Phi vs M Figure Far From GB')
figname=['lowGBD Phi V M Figure flip' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
%% Fit linear equation to Phi vs H -- see lsqfitting for nonlinear
%{

figphiVH=figure;
figphiVH.Name='Phi vs H';
scatter(FPhiGrtest(:)*180/pi,H(:))
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([0.5 6])
title('Declination angle against measured hardness')
figname=['Phi V H Figure Flipped' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')


figphi=figure;
figphi.Name='Phi';
hplot=contourf(X,Y,FPhiGrtest*180/pi(),45,'LineColor','None');
title('Corrected EBSD map')
xlabel('\mum')
ylabel('\mum')
axis image
caxis([0 90])
c=colorbar;
c.Label.String = 'Declination angle /^{o}';
figname=['Phi Figure paper' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')


Hline=H(~isnan(H1));
philine=FPhiGrtest(:)*180/pi;
philine=philine(~isnan(H1));

[hardfit,hardfitq] = polyfit(philine,Hline,1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(hardfit(1)) '*x + ' num2str(hardfit(2))])
hrsq=1 - hardfitq.normr^2 / norm(Hline-mean(Hline))^2;
disp(['R squared is = ' num2str(hrsq)])
% Evaluate fit equation using polyval
[H_est,Hdelta] = polyval(hardfit,philine,hardfitq);
% Add trend line to plot
hold on
plot(philine,H_est,'r','LineWidth',2)
plot(philine,H_est+2*Hdelta,'m--',philine,H_est-2*Hdelta,'m')
scatter(philine,Hline,'bx')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([1.5 5])
%xlim([0 90])
hold off
title('Declination angle against measured hardness with fit')
txt = (['H = ',num2str(hardfit(1),'%4.2f'), '* theta + ', num2str(hardfit(2),'%4.2f')]); %THIS SHIT IS WRONG
text(55, 4.25, txt);%, 'FontSize', 10, 'Color', 'k');
legend({'Linear Fit','90% confidence interval'},'Location','northeast')
figname=['lowGBD Phi V H Figure FlippedFIT' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')



Mline=M(~isnan(M1));
[modfit,modfitq] = polyfit(philine,Mline,1);
% Display evaluated equation y = m*x + b
disp(['Equation is y = ' num2str(modfit(1)) '*x + ' num2str(modfit(2))])
mrsq=1 - modfitq.normr^2 / norm(Mline-mean(Mline))^2;
disp(['R squared is = ' num2str(mrsq)])
% Evaluate fit equation using polyval
[M_est,mdelta] = polyval(modfit,philine,modfitq);
% Add trend line to plot
hold on
plot(philine,M_est,'r--','LineWidth',2)
plot(philine,M_est+2*mdelta,'m--',philine,M_est-2*mdelta,'m--')
scatter(philine,Mline)
xlabel('Phi /degrees')
ylabel('Modulus /GPa')
hold off
title('Phi vs M Figure Far From GB with fit')
figname=['lowGBD Phi V M Figure FlippedFIT' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
%}
%% Picking out individual grains, and looking at where the spread lies. 
%this actually doesn't show much: it's too noisy.
%{

whereg1 = FPhiGrtest > 25*pi/180 & FPhiGrtest < 33*pi/180;
Hg1=H.*whereg1;
Hg1(Hg1==0)=NaN;
FPhiGrtestg1=FPhiGrtest.*whereg1;
FPhiGrtestg1(FPhiGrtestg1==0)=NaN;

%{
scatter(FPhiGrtest(:)*180/pi,Hg1(:),'x','b')
xlabel('Declination Angle /deg')
ylabel('Hardness /GPa')
title('Individual grain Oxygenvs H Figure')
figname=['IndGrain Oxygen V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
%saveas(gcf,fullfile(resultsdir, figname),'png')
%}

figure;
hplot=contourf(X,Y,Hg1,45,'LineColor','None');
title('Isolated EBSD Map')
axis image
c=colorbar;
c.Label.String = 'Oxygen Concentration /arb units';
figname=['FC Figure G4' ebsdname(1:(max(size(ebsdname)-4)))];
%saveas(gcf,fullfile(resultsdir, figname),'png')

%now just look at the bits within that which are high, low, and average
Hg1mean=nanmean(Hg1(:));
Hg1std=nanstd(Hg1(:));

whereg1a = Hg1<Hg1mean-Hg1std;
Hg1a=Hg1.*whereg1a;
Hg1a(Hg1a==0)=NaN;

whereg1b = Hg1>Hg1mean-Hg1std & Hg1<Hg1mean+Hg1std;
Hg1b=Hg1.*whereg1b;
Hg1b(Hg1b==0)=NaN;
whereg1 = Hg1>Hg1mean+Hg1std;
Hg1c=Hg1.*whereg1;
Hg1c(Hg1c==0)=NaN;

hplot=contourf(X,Y,Hg1a,45,'LineColor','None');
%}

%% Low oxygen thresholding, and fitting 
%

whereLowOx = datastack.EPMAO<4;
Hlox=Hmat.*whereLowOx;
whereLowOx = datastack.X<50;
Hlox=Hlox.*whereLowOx;

Hlox(Hlox==0)=NaN;
figphiVHlo=figure;
figphiVHlo.Name='Phi vs H';
scatter(datastack.Phi(:)*180/pi(),Hlox(:))
xlabel('Declination angle  /^{o}')
ylabel('Hardness /GPa')
ylim([2 6])
xlim([0 90])
title({'Declination angle against measured nanoindentation hardness,';' excluding points with high oxygen'})
figname=['Low Oxygen Phi V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

figH=figure;
figH.Name='Hlowox';
contourf(datastack.X,datastack.Y,Hlox,45,'LineColor','None');
title('Nanoindentation map excluding points with high oxygen')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['Hlowox' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
%}




%% Binning the phi data (smearing within grains, spatially blind)

%make an average for all 10deg bins:
phiave=zeros(9,1);
Hphiave=zeros(9,1);
stdevHphiave=zeros(9,1);
countnum=zeros(9,1);
for i = 1:9
    wherei = FPhiGrtest >(i-1)*10*pi/180 & FPhiGrtest <i*10*pi/180;
    phiave(i)=(i-0.5)*10;
    Hphiave(i) = nanmean(H1(wherei));
    stdevHphiave(i) =nanstd(H1(wherei));
    countnum(i)=nnz(wherei);
end
scatter(phiave(:),Hphiave(:))
hold on
errorbar(phiave(:),Hphiave(:),stdevHphiave(:),'LineStyle','none');
xlabel('Phi /degrees')
ylabel('Hardness /GPa')
text(phiave(:)+3,Hphiave(:)+0.5,num2cell(countnum))
title('Low Oxygen Phi vs H Figure Binned to 10deg')
figname=['Low Oxygen Phi V H Figure Binned' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
hold off

%% Pulling out invidual grains of interest

%first lets get rid of bakelite
wheremat = datastack.Y > 642;
wheremat = datastack.X > 170;

wheremat = datastack.H < 12 & datastack.X > 160;
wheremat = datastack.H < 5 & datastack.X > 150;

Hmat=H.*~wheremat;

Hmat=Hmat.*~wheremat;
Hmat(Hmat==0)=NaN;
datastack.Hmat=Hmat;

wheremat = isnan(datastack.Hmat);

datastack.EPMAO=datastack.EPMAO.*~wheremat;
datastack.EPMAO(datastack.EPMAO==0)=NaN;
datastack.EPMABSE=datastack.EPMABSE.*~wheremat;
datastack.EPMABSE(datastack.EPMABSE==0)=NaN;

figure;
hplot=contourf(X,Y,datastack.EPMAOS,45,'LineColor','None');
title('Material only Hardness Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['Hmat Figure paper' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')

whereg1 = datastack.Phi > 55*pi/180 & datastack.Phi < 60*pi/180;
whereg2 = datastack.Phi > 70*pi/180 & datastack.Phi < 75*pi/180;
whereg3 = datastack.Phi > 75*pi/180 & datastack.Phi < 80*pi/180;
whereg4 = datastack.Phi > 80*pi/180 & datastack.Phi < 83*pi/180;
whereg5 = datastack.Phi > 83*pi/180 & datastack.Phi < 86*pi/180;
whereg6 = datastack.Phi > 86*pi/180 & datastack.Phi < 90*pi/180;

Hmat=datastack.Hmat;
Hg1=Hmat.*whereg1;
Hg1(Hg1==0)=NaN;
Hg2=Hmat.*whereg2;
Hg2(Hg2==0)=NaN;
Hg3=Hmat.*whereg3;
Hg3(Hg3==0)=NaN;
Hg4=Hmat.*whereg4;
Hg4(Hg4==0)=NaN;
Hg5=Hmat.*whereg5;
Hg5(Hg5==0)=NaN;
Hg6=Hmat.*whereg6;
Hg6(Hg6==0)=NaN;

Hgtot=cat(3,Hg1, Hg2, Hg3, Hg4, Hg5, Hg6);

figure;
scatter(datastack.EPMAO(:),Hg1(:),'x','b')
hold on
scatter(datastack.EPMAO(:),Hg2(:),'x','k')
scatter(datastack.EPMAO(:),Hg3(:),'x','r')
scatter(datastack.EPMAO(:),Hg4(:),'x','g')
scatter(datastack.EPMAO(:),Hg5(:),'x','m')
scatter(datastack.EPMAO(:),Hg6(:),'x','c')
xlabel('Oxygen conc /arb units')
ylabel('Hardness /GPa')
legend('Phi = 58 deg', 'Phi = 73 deg','Phi = 78 deg','Phi = 81 deg','Phi = 85 deg','Phi = 88 deg');
title('Individual grain Oxygen vs Measured nanoindentation hardness')
figname=['IndGrain Oxygen V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
hold off

for i = 1:6
    Hgi=Hgtot(:,:,i);
    figure;
    hplot=contourf(datastack.X,datastack.Y,Hgi,45,'LineColor','None');
    title('Individual Orientation Hardness Map')
    axis image
    xlabel('\mum')
    ylabel('\mum')
    c=colorbar;
    c.Label.String = 'Hardness /GPa';
    figname=['FC Figure G' num2str(i) ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    delete Hgi
end
%% phi vs GBZ binning



figGBZVH=figure;
figGBZVH.Name='gbz vs H';
scatter(datastack.GBSZ(:),H(:))
xlabel('Declination angle /^{o}')
ylabel('Grain Size /\mum^{2}')
ylim([4 12])
title({'Grain Size against measured hardness'})
figname=['GBZ V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')

%make bins:
binno=15;
gbszave=zeros(binno,1);
Have=zeros(binno,1);
stdevHave=zeros(binno,1);
countnum=zeros(binno,1);
alpha=1; %divvy things up not linearly - but exponent (more in small grains)
beta=log((round(max(max(datastack.GBSZ)))+alpha)/alpha)/binno;
for i = 1:binno
    wherei = datastack.GBSZ > alpha*exp(beta*(i-1)) & datastack.GBSZ <=alpha*exp(beta*i);
    gbszave(i)=alpha*exp(beta*(i-0.5));
    Have(i) = nanmean(H(wherei));
    stdevHave(i) =nanstd(H(wherei));
    countnum(i)=nnz(wherei);
end
scatter(gbszave(:),Have(:))
hold on
errorbar(gbszave(:),Have(:),stdevHave(:),'LineStyle','none');
xlabel('Grain size /\mum^{2}')
ylabel('Hardness /GPa')
text(gbszave(:)+3,Have(:)+0.5,num2cell(countnum))
title('Grain Size vs H Figure Binned')
figname=['Grain Size V H Figure Binned' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')
hold off






%lets try excluding points close to the boundary in big grains? 
%{
binno=10;
gbszave=zeros(binno,1);
Have=zeros(binno,1);
stdevHave=zeros(binno,1);
countnum=zeros(binno,1);
alpha=10; %divvy things up not linearly - but exponent (more in small grains)
beta=log((round(max(max(datastack.GBSZ)))+1)/alpha)/binno;
for i = 1:binno
    wherei = datastack.GBSZ > alpha*exp(beta*(i-1)) & datastack.GBSZ <=alpha*exp(beta*i);
    gbszave(i)=alpha*exp(beta*(i-0.5));
    wherehighGBD = datastack.GBD>sqrt(gbszave(i)/pi*4)*0.1*i/(binno); 
    %also exclude things that are nearer than a certain distance to the GB.
    %estimated distance is by using a circle, halving that radius, and 
    %then loosening the constraint. 
    Hi=datastack.H.*wherehighGBD;
    Hi(Hi==0)=NaN;
    Hii=Hi.*wherei;
    Hii(Hii==0)=NaN;
    Have(i) = nanmean(Hii(:));
    stdevHave(i) =nanstd(Hii(:));
    countnum(i)=nnz(wherei);
    figure; 
    contourf(X,Y,Hii,45,'LineColor','None');
    axis image
end
scatter(gbszave(:),Have(:))
hold on
errorbar(gbszave(:),Have(:),stdevHave(:),'LineStyle','none');
xlabel('Grain size /\mum^{2}')
ylabel('Hardness /GPa')
text(gbszave(:)+3,Have(:)+0.5,num2cell(countnum))
title('Grain Size vs H Figure BinnedCROPPED')
figname=['Grain Size V H Figure BinnedCROPPED' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
hold off
%}
% this isn't so useful - it doesn't quite work i think because it's
% over-sampling in the bulk due to mis-







%% k means?

hplot=contourf(X,Y,Hmat,45,'LineColor','None');

load fisheriris
Xiris = meas(:,3:4);

figure;
plot(Xiris(:,1),Xiris(:,2),'k*','MarkerSize',5);
title 'Fisher''s Iris Data';
xlabel 'Petal Lengths (cm)'; 
ylabel 'Petal Widths (cm)';
rng(1); % For reproducibility
[idx,C] = kmeans(Xiris,3);

%% line profiling

aaatest=improfile(Hmat,[1,1],[1,90])

