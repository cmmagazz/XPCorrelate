%A collection of bits of code to play around with data once you've
%corrected it. I.e. Post-Corrected Data Manipulator.
%CMM 2019

%% GBD Filter 
% Remove points close to a grain boundary from phi vs H plot to see the
% effect of thsi
wherehighGBD = datastack.GBD>7;
Hhgbd=datastack.H.*wherehighGBD;
Hhgbd(Hhgbd==0)=NaN;
wherelowGBD = datastack.GBD<2;
Hlgbd=datastack.H.*wherelowGBD;
Hlgbd(Hlgbd==0)=NaN;

%number of points
sum(wherehighGBD(:))

figphiVHlgbd=figure;
figphiVHlgbd.Name='Phi vs H';
scatter(datastack.Phirefl(:)*180/pi(),Hhgbd(:) ,'x')
xlabel('Declination angle  /^{o}')
ylabel('Hardness /GPa')
ylim([1 5])
xlim([0 90])
title({'Declination angle against measured nanoindentation hardness,';' excluding points near grain boundaries'})
figname=['LowGBD Phi V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

XPCcontourf(Hhgbd,'title',"Nanoindentation hardness - excluding points near GB",...
    'cunits',"Hardness /GPa",'saveq',1)


figphiVHlgbd=figure;
figphiVHlgbd.Name='Phi vs H';
scatter(datastack.Phirefl(:)*180/pi,Hhgbd(:))
hold on
scatter(datastack.Phirefl(:)*180/pi,Hlgbd(:),'rx')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
legend({'Points higher than 7um','Points closer than 2um'})
ylim([1.5 5])
title({'Declination angle against measured hardness';'displaying only points close to grain boundaries'})
figname=['lowGBD Phi V H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')
XPCcontourf(Hhgbd)
XPCcontourf(Hlgbd)

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
contourf(X,Y,Hhgbd,45,'LineColor','None');
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

%% Picking out individual grains/angles, and looking at where the spread lies. 
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

%% Low oxygen thresholding 

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

%% Binning the phi data (smearing within grains, spatially blind)
%make an average for all 10deg bins:
phiave=zeros(9,1);
Hphiave=zeros(9,1);
stdevHphiave=zeros(9,1);
countnum=zeros(9,1);
for i = 1:9
    wherei = FPhiGrtest >(i-1)*10*pi/180 & FPhiGrtest <i*10*pi/180;
    phiave(i)=(i-0.5)*10;
    Hphiave(i) = nanmean(Hhgbd(wherei));
    stdevHphiave(i) =nanstd(Hhgbd(wherei));
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

%% Pulling out invidual areas of interest
%Random bits of code to just pull out specific areas
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


%% Declination angle vs GB area binning

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
