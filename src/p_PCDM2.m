%% PCDM2
%Script with even more PCDM, largely following response to Angus'
%comments from the first dialogue of manuscripts (Mar2020)

%% Phi vs H with colour based on GBD - original dataset

figure;
%scatter(datastack.Phi(:)*180/pi(),datastack.H(:),pointsize,datastack.GBD(:),'filled');
temp=[datastack.Phireflseccorr(:)*180/pi(),datastack.H(:),datastack.GBDseccorr(:),10*log(datastack.GBDseccorr(:))];
temp=sortrows(temp,3);
scatter(temp(:,1),temp(:,2),[],temp(:,3),'filled');
c=colorbar();
c.Label.String = 'Distance to GB /\mum';
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5*nanstd(H1(:))])
%xlim([0 90])
title({'Declination angle against measured hardness',' with distance to grain boundary'})
figname=['Phi V H Figure_colours' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end
%% Normalising H and plotting w.r.t GBD
Hpred=a(1)*cos(2*datastack.Phireflseccorr)+a(2);

datastack.Hsmoothed=smoothdata(datastack.H,1,'gaussian',2);
datastack.Hsmoothed=smoothdata(datastack.Hsmoothed,2,'gaussian',2);

Hnormalisedphi=datastack.Hsmoothed-Hpred;

XPCcontourf(Hpred,'title',"Expected Hardness",'cunits',"Hardness /GPa",'saveq',1)
XPCcontourf(Hnormalisedphi,'title',"Normalised Hardness",'cunits',"Hardness /GPa",'saveq',1)

%the famous trumpet graph
figure;
scatter(datastack.GBDseccorr(:),Hnormalisedphi(:));
xlabel('Distance to GB /\mum')
ylabel('Normalised Hardness /GPa')
%ylim([nanmean(H1(:))-5*nanstd(H1(:)) nanmean(H1(:))+5*nanstd(H1(:))])
%xlim([0 90])
ylim([-2 2])
title({'GB distance against normalised measured hardness'})
figname=['GBD V H norm figure blankseccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end


%% Binning the GBD and averaging. 
%make an average for all 10deg bins:
hgbdbinno=50;
gbdave=zeros(hgbdbinno,1);
Hgbdave=zeros(hgbdbinno,1);
stdevHgbdave=zeros(hgbdbinno,1);
countnum=zeros(hgbdbinno,1);
for i = 1:hgbdbinno
    wherei = datastack.GBDseccorr >(i-1)*max(datastack.GBDseccorr(:))/hgbdbinno & datastack.GBDseccorr <i*max(datastack.GBDseccorr(:))/hgbdbinno;
    gbdave(i)=(i-0.5)*max(datastack.GBDseccorr(:))/hgbdbinno;
    Hgbdave(i) = nanmean(Hnormalisedphi(wherei));
    stdevHgbdave(i) =nanstd(Hnormalisedphi(wherei));
    countnum(i)=nnz(wherei);
end
scatter(gbdave(:),Hgbdave(:))
hold on
errorbar(gbdave(:),Hgbdave(:),stdevHgbdave(:),'LineStyle','none');
xlabel('Distance to GB /\mum')
ylabel('Hardness /GPa')
%text(gbdave(:)+3,Hgbdave(:)+0.05,num2cell(countnum))
title('Hardness near GB, Averaged')
figname=['HnearGB Binned seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
hold off


%% Visualising the difference between Phi predicted and Phi, first and second correction
% Another way to visualise this: 
imgcomp=imshowpair(rot90(datastack.Phipred),rot90(datastack.Phi))
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))

imgcomp=imshowpair(rot90(datastack.Phipred),rot90(datastack.Phi),'diff')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi2' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))

%Now doing this with the second correction:
imgcomp=imshowpair(rot90(datastack.Phipredseccorr),rot90(datastack.Phiseccorr),'diff')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
figname=['Phipred vs Phi2seccorr' ebsdname(1:(max(size(ebsdname)-4))) '.png'];
imwrite(imgcomp.CData,fullfile(resultsdir, figname))

