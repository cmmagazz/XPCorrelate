function f_plotfig(datastack,resultsdir,ebsdname,saveasfigq)

%Plot fig - CMM script to plot out the various graphs from the data. 

resolution=['-r' num2str(600)];

meanH=nanmean(datastack.H(:));
stdH=nanstd(datastack.H(:));
meanM=nanmean(datastack.M(:));
stdM=nanstd(datastack.M(:));

figure;
try
    meanH=nanmean(datastack.Hmat(:));
    stdH=nanstd(datastack.Hmat(:));
    hplot=contourf(datastack.X,datastack.Y,datastack.Hmat,455,'LineColor','None');
catch
    hplot=contourf(datastack.X,datastack.Y,datastack.H,455,'LineColor','None');
end
if meanH>stdH
    caxis([meanH-2*stdH meanH+2*stdH])
else 
    caxis([min(hplot(:)) meanH+1*stdH])
end
title('Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness (GPa)';
figname=['Hardness Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end

% smoothing currently not implemented, but here if need be.
%H2=smoothdata(H,2,'gaussian',7);
%H2=smoothdata(H2,1,'gaussian',7);
%{
figure;
hplot=contourf(datastack.X,datastack.Y,datastack.M,455,'LineColor','None');
title('Modulus')
xlabel('\mum')
ylabel('\mum')
axis image
caxis([meanM-4*stdM meanM+4*stdM])
c=colorbar;
c.Label.String = 'Modulus (GPa)';
figname=['Modulus Figure ' filename(1:(max(size(filename)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 saveas(gcf,fullfile(resultsdir, figname),'fig') end
%}

figphi=figure;
figphi.Name='Phi';
hplot=contourf(datastack.X,datastack.Y,datastack.Phirefl*180/pi(),45,'LineColor','None');
title('Registered EBSD map')
xlabel('\mum')
ylabel('\mum')
axis image
zlim([0 90])
c=colorbar;
c.Label.String = 'Declination angle /^{o}';
figname=['Phi Figure paper' ebsdname(1:(max(size(ebsdname)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end


figbc=figure;
figbc.Name='BC';
hplot=contourf(datastack.X,datastack.Y,datastack.BCebsd,45,'LineColor','None');
title('Registered EBSD map - BC')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'BC /arb units';
figname=['BC Figure paper' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end


figphiVH=figure;
figphiVH.Name='Phi vs H';
scatter(datastack.Phirefl(:)*180/pi,datastack.H(:),'x')
title('Declination angle against measured nanoindentation hardness')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(datastack.H(:))-5*nanstd(datastack.H(:)) nanmean(datastack.H(:))+5*nanstd(datastack.H(:))])
xlim([0 90])
figname=['Phi vs H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
print(fullfile(resultsdir, figname),'-dpng',resolution)
if saveasfigq==1 
    saveas(gcf,fullfile(resultsdir, figname),'fig') 
end

%{
figphiVM=figure;
figphiVM.Name='Phi vs M';
scatter(datastack.Phi(:)*180/pi,M(:))
xlabel('Phi /degrees')
ylabel('Modulus /GPa')
figname=['Phi vs M Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
if saveasfigq==1 saveas(gcf,fullfile(resultsdir, figname),'fig') end
%}

close all 

end
