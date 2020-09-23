function [datastack] = f_fixEBSDdistortion(ebsd,X,Y,fullres,microscope,primphase, EBSDREFq,resultsdir)
% Script for correcting the distortion in the EBSD map based on what we
% presume to be a reliable hardness map. Selection of points are made 
% (at least 4) in both maps, and an affine transformation that links the
% hardness points onto the ebsd map is found. All coordinates on the 
% hardness map are then transformed to form "query coordinates"
% on the ebsd map. The result is then gridded, and outputted to datastack
% which contains the X, Y, Express, and EBSD data in the same size arrays.

%CMM and IH and AJW 2020.

%Datastack structure: 
% datastack.X               %X position
% datastack.Y               %Y position
% datastack.S               %Surface displacement
% datastack.D               %Depth
% datastack.L               %Load
% datastack.M               %Modulus
% datastack.St              %Stiffness^2/Load
% datastack.H               %Hardness
% datastack.phi1            %phi1  (ebsd)
% datastack.Phi             %Phi   (ebsd)
% datastack.phi2            %phi2  (ebsd)
% datastack.phase           %phase (ebsd)

%% hardness input
%First make sure the data is appropriately structured. FROM HERE ON OUT,
%data should be structured in grids which start at 0,0, and X(2,1) is
%larger than X(1,1), and Y(1,2) is larger than Y(1,1). This goes against
%Matlab's conventional Row Column notation, but follows the cartesian form
%to relate position in matrices to positions on the cartesian grid
%described by X and Y. 

%quick check:
i=1;
while isnan(X(i,1)) || isnan(X(i+1,1))
    i=i+1;
end
if X(i+1,1)<X(i,1)
    X=flipud(X);
    Y=flipud(Y);
    fullres=flipud(fullres);
end
i=1;
while isnan(Y(1,i)) || isnan(Y(1,i+1))
    i=i+1;
end
if Y(1,i+1)<Y(1,i)
    X=fliplr(X);
    Y=fliplr(Y);
    fullres=fliplr(fullres);
end

if X(1,1)==0 && Y(1,1)==0 
    disp('H zerod correctly')
else
    X= X-min(min(X));
    Y= Y-min(min(Y));
end

% Himage
H=fullres(:,:,6);%Extract hardness
H(H>1000)=0;
figure
hplot=contourf(X,Y,H,45,'LineColor','None');
axis image
caxis([nanmean(H(:))-1*nanstd(H(:)) nanmean(H(:))+1*nanstd(H(:))])
title('Reference Selection in Hardness(>4)')
[x_og,y_og] = getpts; %obtain the "og" coordinates: the original, absolute reference frame for all future points.
close all

%plot it again with the points on for helpful guide, and save.
figure
hplot=contourf(X,Y,H,45,'LineColor','None');
hold on
scatter(x_og,y_og,'rx')
for i=1:size(x_og,1)
    text(x_og(i)+0.01*max(X(:)), y_og(i)+0.01*max(Y(:)),num2str(i));   %put up the text at the given locations
end
title('References Selected On Nanoindentation Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness /GPa';
figname=['refselect_H figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')
%% ebsd input 
% EBSD settings and cleaning:
%

if strcmp(primphase,'default')==1
    primphase=char(ebsd.mineralList(2));
end
setMTEXpref('xAxisDirection','west');

if EBSDREFq==0 || EBSDREFq==2
    ebsd=ebsd(primphase); 
end
%clean things up
[Grains,ebsd.grainId] = calcGrains(ebsd,'boundary','convexhull','angle',10*degree); %calc grains - this is useful for cleaning up
large_grains = Grains(Grains.grainSize >= 3); 
ebsd_clean = ebsd(large_grains);
%fill EBSD
% delete nonindexed and fill them with nn existing phase 
ebsd_clean('n')=[] ;
ebsd_clean = fill(ebsd_clean) ;
ebsd=ebsd_clean;
%% EBSD manipulation
%We deal with EBSD data in the same way as hardness: based on 0,0 and going
%up in the correct way.
[ebsdGrid] = gridify(ebsd);
% now set the x and y coordinates as what they are in the ebsd map
% (cropped and cleaned)
x_ebsd = ebsdGrid.prop.x;
y_ebsd = ebsdGrid.prop.y;
try
    BCebsd = ebsdGrid.prop.bc;
catch
    BCebsd = ebsdGrid.prop.RadonQuality;
end

        
%extract the things we care about: angles
if EBSDREFq==0 || EBSDREFq==2 
    phi1 = ebsdGrid.orientations.phi1;
    Phi  = ebsdGrid.orientations.Phi;
    phi2 = ebsdGrid.orientations.phi2;
elseif EBSDREFq==1
    allphi1s=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);
    allPhis=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);
    allphi2s=zeros(size(x_ebsd,1),size(y_ebsd,2),length(ebsdGrid.mineralList)-1);

    for i = 1:length(ebsdGrid.mineralList)-1 %go through all the phases and get their orientations
        try
            phasenebsd=gridify(ebsdGrid(ebsdGrid.mineralList(i+1)));
        catch
            i=i+1;
            if i<=length(ebsdGrid.mineralList)-1
                phasenebsd=gridify(ebsdGrid(ebsdGrid.mineralList(i+1)));
            end
        end
        allphi1s(:,:,i) = phasenebsd.orientations.phi1;
        allPhis(:,:,i)  = phasenebsd.orientations.Phi;
        allphi2s(:,:,i) = phasenebsd.orientations.phi2;
    end
    phi1=mean(allphi1s,3,'omitnan'); delete allphi1s %get the phi of the relevant phase
    Phi=mean(allPhis,3,'omitnan'); delete allPhis %don't worry about mean
    phi2=mean(allphi2s,3,'omitnan'); delete allphi2s %if it's not that phase it's nan and this ignores nan
end
phase= ebsdGrid.phase;
%check: are the x and y values changing in the same way as the
%hardness? i.e. does x(2,1) correspond to the pixel on the right of x(1,1)
if x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1) %if everythings fine
    %do nothing
elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
    x_ebsd=flipud(x_ebsd);
    y_ebsd=flipud(y_ebsd);
    phi1=flipud(phi1);
    Phi=flipud(Phi);
    phi2=flipud(phi2);
    phase=flipud(phase);
    BCebsd=flipud(BCebsd);
elseif x_ebsd(2,1)> x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
    x_ebsd=fliplr(x_ebsd);
    y_ebsd=fliplr(y_ebsd);
    phi1=fliplr(phi1);
    Phi=fliplr(Phi);
    phi2=fliplr(phi2);
    phase=fliplr(phase);
    BCebsd=fliplr(BCebsd);
elseif x_ebsd(2,1)<= x_ebsd(1,1) && y_ebsd(1,2)<=y_ebsd(1,1)
    x_ebsd=transpose(x_ebsd);
    y_ebsd=transpose(y_ebsd);
    phi1=transpose(phi1);
    Phi=transpose(Phi);
    phi2=transpose(phi2);
    phase=transpose(phase);
    BCebsd=transpose(BCebsd);
end

if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 
    %do nothing
else
    x_ebsd= x_ebsd-min(min(x_ebsd));
    y_ebsd= y_ebsd-min(min(y_ebsd));
end

%HOWEVER, each ebsd system will decide which was this is. I've tried to get
%the right info for each microscope, but IF THIS ISN'T RIGHT JUST CLOSE THE
%REF FIGURE AND IT WILL MANUALLY LET YOU PICK THE DIRECTION
% put the first point on (0, 0)
if strcmp(microscope,'merlin') 
    %if it's a merlin:
    phi1=flipud(fliplr(phi1));
    Phi=flipud(fliplr(Phi));
    phi2=flipud(fliplr(phi2));
    phase=flipud(fliplr(phase));
    BCebsd=flipud(fliplr(BCebsd));
elseif strcmp(microscope,'evo') 
    %if it's a evo:
    phi1=fliplr(phi1);
    Phi=fliplr(Phi);
    phi2=fliplr(phi2);
    phase=fliplr(phase);
    BCebsd=fliplr(BCebsd);
elseif strcmp(microscope,'tescan') 
    %if it's oxford instruments
    phi1=fliplr(phi1);
    Phi=fliplr(Phi);
    phi2=fliplr(phi2);
    phase=fliplr(phase);
    BCebsd=fliplr(BCebsd);
elseif strcmp(microscope,'xbeam') 
    %if it's the crossbeam instruments
    %the xbeam is frustrating
    [x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd]=fixXbeam(x_ebsd,y_ebsd,phi1, phi2, Phi, phase,BCebsd);
end

%quick check:
if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 && x_ebsd(2,1)>x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
    disp('ebsd zerod correctly')
else
    error('ebsd data not correctly structured')
end
%
%% select points on the manipulated ebsd
% IF THIS IS THE WRONG WAY, just close the figure
try
    ebsdfig=figure;
    if EBSDREFq==0 
        hplot=contourf(x_ebsd,y_ebsd,Phi,45,'LineColor','None');
    elseif EBSDREFq==1
        hplot=contourf(x_ebsd,y_ebsd,phase,45,'LineColor','None');
    elseif EBSDREFq==2
        hplot=contourf(x_ebsd,y_ebsd,BCebsd,45,'LineColor','None');
        caxis([nanmean(BCebsd(:))-nanstd(BCebsd(:)) nanmean(BCebsd(:))+nanstd(BCebsd(:))]) 
    end
    hold on
    axis image
    hold off
    title('Reference Selection in EBSD (>4)')
    [x_transREF,y_transREF] = getpts;
    close(ebsdfig)
catch
    [x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,x_transREF,y_transREF]=f_ebsdwrongway(x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,EBSDREFq);
end

%plot it again with the points on for helpful guide, and save.
figure
if EBSDREFq==0 
    hplot=contourf(x_ebsd,y_ebsd,Phi,45,'LineColor','None');
    c=colorbar;
    c.Label.String = 'Declination angle /^{o}';
elseif EBSDREFq==1
    hplot=contourf(x_ebsd,y_ebsd,phase,45,'LineColor','None');
    c=colorbar;
    c.Label.String = 'Phase';
elseif EBSDREFq==2
    hplot=contourf(x_ebsd,y_ebsd,BCebsd,45,'LineColor','None');
    caxis([nanmean(BCebsd(:))-nanstd(BCebsd(:)) nanmean(BCebsd(:))+nanstd(BCebsd(:))]) 
    c=colorbar;
    c.Label.String = 'Band Contrast /arb units';
end
hold on
scatter(x_transREF,y_transREF)
for i=1:size(x_transREF,1)
    text(x_transREF(i)+0.01*max(X(:)), y_transREF(i)+0.01*max(Y(:)),num2str(i));   %put up the text at the given locations
end
title('References Selected On EBSD Map')
xlabel('\mum')
ylabel('\mum')
axis image
figname=['refselect_EBSD figure'];
saveas(gcf,fullfile(resultsdir, figname),'png')
close all
%% transformation calculation

movingPoints  = [x_og    y_og   ];
fixedPoints = [x_transREF y_transREF];

transtype='affine'; %select what type of transformation you want
tform = fitgeotrans(movingPoints,fixedPoints,transtype);
% transform the data
A = tform.T;

%transform the x and y coordinates of the nanoindentation data
%using the affine transformation
[locebsdcorr(:,1),locebsdcorr(:,2)] = transformPointsForward(tform,X(:),Y(:));

%% interpolation section
%using griddedInterpolant, setup a variable that contains all the
%information, and then probe at the locations of the distorted x and y
%locations from the nanoindentation data
% 'nearest' is used as we don't want to invent orientations

 phi1interp = griddedInterpolant(x_ebsd,y_ebsd,phi1,'nearest');
 phi1new = phi1interp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 Phiinterp = griddedInterpolant(x_ebsd,y_ebsd,Phi,'nearest');
 Phinew = Phiinterp(locebsdcorr(:,1), locebsdcorr(:,2));
  
 phi2interp = griddedInterpolant(x_ebsd,y_ebsd,phi2,'nearest');
 phi2new = phi2interp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 phaseinterp = griddedInterpolant(x_ebsd,y_ebsd,phase,'nearest');
 phasenew = phaseinterp(locebsdcorr(:,1), locebsdcorr(:,2));
 
 BCebsdinterp = griddedInterpolant(x_ebsd,y_ebsd,BCebsd,'nearest');
 BCebsdnew = BCebsdinterp(locebsdcorr(:,1), locebsdcorr(:,2));
 
%Gridify into a matrix 
phi1newG = f_gridify_vector(phi1new,size(X,1),size(Y,2))';
PhinewG = f_gridify_vector(Phinew,size(X,1),size(Y,2))';
phi2newG = f_gridify_vector(phi2new,size(X,1),size(Y,2))';
phasenewG = f_gridify_vector(phasenew,size(X,1),size(Y,2))';
BCebsdnewG = f_gridify_vector(BCebsdnew,size(X,1),size(Y,2))';

%output into a class datastack with reasonable names
datastack.X     = X;                 %X position
datastack.Y     = Y;                 %Y position
datastack.S     = fullres(:,:,1);   %Surface displacement
datastack.D     = fullres(:,:,2);   %Depth
datastack.L     = fullres(:,:,3);   %Load
datastack.M     = fullres(:,:,4);   %Modulus
datastack.St    = fullres(:,:,5);   %Stiffness^2/Load
datastack.H     = fullres(:,:,6);   %Hardness
datastack.phi1  = phi1newG;          %phi1  (ebsd)
datastack.Phi   = PhinewG;           %Phi   (ebsd)
datastack.phi2  = phi2newG;          %phi2  (ebsd)
datastack.phase = phasenewG;         %phase (ebsd)
datastack.BCebsd= BCebsdnewG;        %Band Contrast (ebsd)

end
