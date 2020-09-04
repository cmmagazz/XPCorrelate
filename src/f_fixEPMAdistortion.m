function [datastack]=f_fixEPMAdistortion(datastack,C,XC,YC,cropq, EPMArefq)

% Script for correcting the distortion in the EPMA map based on a reliable
% hardness map. See fix_EBSDdistortion for details
% Output is fullstack, which contains the X,Y,Express,EBSD, and EPMA all in
% the same size and related. CMM 2019.



%% hardness input
%First make sure the data is appropriately structured. FROM HERE ON OUT,
%data should be structured in grids which start at 0,0, and X(2,1) is
%larger than X(1,1), and Y(1,2) is larger than Y(1,1). This goes against
%Matlab's conventional Row Column notation, but follows the cartesian form
%to relate position in matrices to positions on the cartesian grid
%described by X and Y. 

%quick check:

% Himage
H=datastack.H;
X=datastack.X;
Y=datastack.Y;
figure
hplot=contourf(X,Y,H,455,'LineColor','None');
axis image
caxis([nanmean(H(:))-1*nanstd(H(:)) nanmean(H(:))+1*nanstd(H(:))])
title('Reference Selection in Hardness(>4)')
[x_og,y_og] = getpts; %obtain the "og" coordinates: the original, absolute reference frame for all future points.

%% cropping of EPMA data
if cropq==1

    figure
    hplot=contourf(XC,YC,C,45,'LineColor','None');
    title('Cropping EPMA map (select corners)')
    axis image
    [xc_crop,yc_crop] = getpts;
    close all

    C=C(round(min(xc_crop)):round(max(xc_crop)),round(min(yc_crop)):round(max(yc_crop)));
    XC=XC(round(min(xc_crop)):round(max(xc_crop)),round(min(yc_crop)):round(max(yc_crop)));
    YC=YC(round(min(xc_crop)):round(max(xc_crop)),round(min(yc_crop)):round(max(yc_crop)));
end


%% section to align to a point on EPMA, and then interp using the hardness grid
%start at 0,0
XC=XC-min(min(XC));
YC=YC-min(min(YC));

figure
if EPMArefq==0 && size(C,3)==2
    hplot=contourf(XC,YC,C(:,:,1),45,'LineColor','None');
elseif EPMArefq==1 && size(C,3)==2
    hplot=contourf(XC,YC,C(:,:,2),45,'LineColor','None');
else % if there isn't a BSE image, then just use the standard irregardless of EPMArefq
    hplot=contourf(XC,YC,C,45,'LineColor','None');
end 
%{
if nanmean(hplot(:))>nanstd(hplot(:))
    caxis([nanmean(hplot(:))-0.01*nanstd(hplot(:)) nanmean(hplot(:))+0.01*nanstd(hplot(:))])
else 
    caxis([min(hplot(:)) nanmean(hplot(:))+nanstd(hplot(:))])
end
%} 
%fix this
caxis([95 105])
hold on
axis image
hold off
title('Reference Selection in EPMA (>4)')
[x_transREF,y_transREF] = getpts;
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
[locepmacorr(:,1),locepmacorr(:,2)] = transformPointsForward(tform,X(:),Y(:));

%% interpolation section
%using griddedInterpolant, setup a variable that contains all the
%information, and then probe at the locations of the distorted x and y
%locations from the nanoindentation data

if size(C,3)==1 %if there's no BSE image
    Cinterp = griddedInterpolant(XC,YC,C,'linear');
    Onew = Cinterp(locepmacorr(:,1), locepmacorr(:,2));
    OnewG = gridify_vector(Onew,size(X,1),size(Y,2))';%Gridify into a matrix 
    
elseif size(C,3)==2 % if there is
    Cinterp = griddedInterpolant(XC,YC,C(:,:,1),'linear');
    Onew = Cinterp(locepmacorr(:,1), locepmacorr(:,2));
    OnewG = gridify_vector(Onew,size(X,1),size(Y,2))';%Gridify into a matrix 

    BSEinterp = griddedInterpolant(XC,YC,C(:,:,2),'linear');
    BSEnew = BSEinterp(locepmacorr(:,1), locepmacorr(:,2));
    BSEnewG = gridify_vector(BSEnew,size(X,1),size(Y,2))';%Gridify into a matrix 

    datastack.EPMAO     = OnewG;    %EPMA map O
    datastack.EPMABSE   = BSEnewG;  %BSE frmo the EPMA
end
end