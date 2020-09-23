function [datastack,ebsdCorrectedLGG,Grains]=g_dist2grainb_secondcorr(resultsdir, ebsdname,datastack,currdate,ebsd)
%script to calculate, at every point, the distance to the grain boundary.
%inputs are just the directory and filename of the CORRECTED and 
%appropriately sized ebsd map. import into mtex, and then distances 
%are calculated thereafter. 

try
    CS=f_phaseTranslate(ebsd);
catch
    CS=loadCIF('Ti-Titanium-alpha');
end

ebsdCorrected = loadEBSD(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'secondcorrected' currdate]), CS,...
'interface','generic','ColumnNames', { 'x' 'y' 'phase' 'phi1' 'Phi' 'phi2'}, 'Columns', [1 2 3 4 5 6], 'Bunge');
primphase=char(ebsdCorrected.mineralList(2));

%ebsdCorrected = rot90(ebsdCorrected,1);
setMTEXpref('xAxisDirection','east');%FLIPPING HERE - FIX THIS IT@S NONSENSE
%clean things up
ebsdCorrected('n')=[] ;
[Grains,ebsdCorrected.grainId] = calcGrains(ebsdCorrected,'boundary','convexhull','angle',5*degree); %calc grains - this is useful for cleaning up
%{
%large_grains = Grains(Grains.grainSize >= 20); 
%ebsdCorrectedLG = ebsdCorrected(large_grains);
%[Grains,ebsdCorrectedLG.grainId] = calcGrains(ebsdCorrectedLG,'boundary','convexhull','angle',5*degree);
%
%}
ebsdCorrected = fill(ebsdCorrected) ;
ebsdCorrected = smooth(ebsdCorrected('indexed'),'fill',Grains);
ebsdCorrectedLGG=gridify(ebsdCorrected); 
[Grains,ebsdCorrectedLGG.grainId] = calcGrains(ebsdCorrectedLGG,'boundary','convexhull','angle',5*degree); %calc grains - this is useful for cleaning up
Grains = smooth(Grains,5);


plot(ebsdCorrectedLGG(primphase),ebsdCorrectedLGG(primphase).orientations)
figname=['EBSDCorrectedIPF_seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
close all
%{
ebsdt=flipud(ebsd);
plot(ebsdCorrected(primphase),ebsdCorrected(primphase).orientations)
plot(ebsd(primphase),ebsd(primphase).orientations)
figname=['EBSDOriginalIPF' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
%}

exportCTF(ebsdCorrectedLGG,fullfile(resultsdir,'ebsdseccorr.ctf'))


GBD=zeros(size(ebsdCorrectedLGG.prop.x,1),size(ebsdCorrectedLGG.prop.y,2),1);
GBSZ=zeros(size(ebsdCorrectedLGG.prop.x,1),size(ebsdCorrectedLGG.prop.y,2),1);
gID=zeros(size(ebsdCorrectedLGG.prop.x,1),size(ebsdCorrectedLGG.prop.y,2),1);
gmo=zeros(size(ebsdCorrectedLGG.prop.x,1),size(ebsdCorrectedLGG.prop.y,2),3);

for i=1:size(ebsdCorrectedLGG.prop.x,1)
    for j=1:size(ebsdCorrectedLGG.prop.y,2)
        g=ebsdCorrectedLGG.prop.grainId(i,j);
        %just need to look at every point in ebsdCorrectedLGG and it's
        %respective grainID, and compare THAT to the closest point in the
        %grains(i).boundary.
        gbI=Grains(g).boundary('indexed');
        distances = sqrt((gbI.x-ebsdCorrectedLGG.prop.x(i,j)).^2 + (gbI.y-ebsdCorrectedLGG.prop.y(i,j)).^2);
        [shortd,locshortd] = min(distances);
        GBD(i,j)=shortd;
        GBSZ(i,j)=Grains(g).area; %measured in square microns
        gID(i,j)=g;
        %get the misorientations too 
        try
            gmo(i,j,:)=[gbI(locshortd).misorientation.phi1 gbI(locshortd).misorientation.Phi gbI(locshortd).misorientation.phi2];
        catch
            gmo(i,j,:)=[0 0 0];
        end
    end
    f=waitbar(i / size(ebsdCorrectedLGG.prop.x,1));
end
close(f)

X=datastack.X;
Y=datastack.Y;

if size(GBD,2)==size(X,1) && size(GBD,1)==size(X,2)
    GBD=GBD';
    GBSZ=GBSZ';
    gID=gID';
    gmo=permute(gmo,[2 1 3]);
end

if size(GBD)==size(X) %check if dimensions agree, else interpolate using X,Y loc
else
    ebsdlocx = f_gridify_vector(ebsdCorrectedLGG.prop.x,size(GBD,1),size(GBD,2))';
    ebsdlocy = f_gridify_vector(ebsdCorrectedLGG.prop.y,size(GBD,1),size(GBD,2))';

    GBDinterp = griddedInterpolant(ebsdlocx',ebsdlocy',GBD','nearest');
    GBDnew = GBDinterp(-X(:), -Y(:));
    GBD = f_gridify_vector(GBDnew,size(X,1),size(Y,2))';

    GBSZinterp = griddedInterpolant(ebsdlocx',ebsdlocy',GBSZ','nearest');
    GBSZnew = GBSZinterp(-X(:), -Y(:));
    GBSZ = f_gridify_vector(GBSZnew,size(X,1),size(Y,2))';

    gIDinterp = griddedInterpolant(ebsdlocx',ebsdlocy',gID','nearest');
    gIDnew = gIDinterp(-X(:), -Y(:));
    gID = f_gridify_vector(gIDnew,size(X,1),size(Y,2))';
    
    gmointerp = griddedInterpolant(ebsdlocx',ebsdlocy',gmo','nearest');
    gmonew = gmointerp(-X(:), -Y(:));
    gmo = f_gridify_vector(gmonew,size(X,1),size(Y,2))';
end


datastack.GBDseccorr=GBD;
datastack.GBSZseccorr=GBSZ;
datastack.gIDseccorr=gID;
datastack.gmoseccorr=gmo;

XPCcontourf(sqrt(gmo(:,:,1).^2+gmo(:,:,2).^2+gmo(:,:,3).^2)./GBD)
%
figGBD=figure;
figGBD.Name='Grain Boundary Distances - Second Correction';
hplot=contourf(datastack.X,datastack.Y,datastack.GBDseccorr,45,'LineColor','None');
title(figGBD.Name)
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Distance to GB /\mum';
figname=['Grain Boundary Distances Figure SECCORR' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

figGBsz=figure;
figGBsz.Name='Grain Boundary Area - Second Correction';
hplot=contourf(datastack.X,datastack.Y,datastack.GBSZseccorr,45,'LineColor','None');
title(figGBsz.Name)
xlabel('\mum')
ylabel('\mum')
c=colorbar;
c.Label.String = 'Grain Area /\mum^{2}';
axis image
figname=['Grain Boundary Area Figure seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

figGID=figure;
figGID.Name='Grain ID - Second Correction';
hplot=contourf(datastack.X,datastack.Y,datastack.gIDseccorr,max(max(datastack.gIDseccorr)));
title(figGID.Name)
c=colorbar;
c.Label.String = 'Grain ID';
xlabel('\mum')
ylabel('\mum')
axis image
figname=['Grain ID Figure seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

H=datastack.H;
H(H>1000)=NaN;
figHVPHIVGBD=figure;
figHVPHIVGBD.Name='GBD vs phi vs H';
try
    scatter3(datastack.GBDseccorr(:),datastack.Phiseccorr(:)*180/pi,datastack.Hmat(:))
catch
    scatter3(datastack.GBDseccorr(:),datastack.Phiseccorr(:)*180/pi,datastack.H(:))
end
xlabel('GBD')
ylabel('Phi /degrees')
zlabel('Hardness /GPa')
figname=['GBD vs phi vs H Figure seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')


figHVPHIVGBSZ=figure;
figHVPHIVGBSZ.Name='GBSZ vs phi vs H';
try
    scatter3(datastack.GBSZseccorr(:),datastack.Phiseccorr(:)*180/pi,datastack.Hmat(:))
catch
    scatter3(datastack.GBSZseccorr(:),datastack.Phiseccorr(:)*180/pi,datastack.H(:))
end
xlabel('GBSZ')
ylabel('Phi /degrees')
zlabel('Hardness /GPa')
figname=['GBSZ vs phi vs H Figure seccorr' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

close all
%}
end
