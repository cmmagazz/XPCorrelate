function g_secondcorrection(datastack,buffersize)
% This function performs a second correction to the Phi values. The way it
% works is by assuming that errors arise from grain boundaries being
% misplaced. The values near the grain boundaries are filtered, with
% distance set by BUFFERSIZE, and a fitting similar to lsqfitting is
% performed. The fit is then used to "predict" a phi for every point in the
% hardness map. This should not be nominally correct, but should have the
% grain boundaries in the right place. 

% We then edge filter it, and will end up with two plots: the predicted phi
% edge map, and the first corrected ebsd phi edge map. These should "look"
% the same. I.e. their bright areas should be in the same place, as
% boundaries should lie on boundaries. These two images are fed into a
% an ECC script from an external source, which should be added to the path
% here: 
% https://www.mathworks.com/matlabcentral/fileexchange/27253-ecc-image-alignment-algorithm-image-registration
% I've put this in the external folder. 
addpath(genpath('external'));
% Their guide is much better to what mine could be in how it works. 
% We then save this as an ebsd, re-run the grain boundary distance
% analysis, and spit out some useful graphs. This also spits out the
% following: 
% 
% datastack.Phiseccorr             %Phi, second corrected
% datastack.phi1seccorr            %phi1, second corrected
% datastack.phi2seccorr            %phi2, second corrected
% datastack.phaseseccorr           %phase, second corrected
% datastack.BCeBCebsdseccorrbsd    %Band Contrast, second corrected

if exist('buffersize','var')~=1
    buffersize=input('What buffer size should be used?');
end

%first do a fitting to get predicted values
wherehighGBD = datastack.GBD>buffersize;
H1=datastack.H.*wherehighGBD;
H1(H1==0)=NaN;
%turn into column of data and remove nan
xdata = datastack.Phirefl(:)*180/pi;
H2fit=H1(:);
xdata=xdata(~isnan(H2fit));
H2fit=H2fit(~isnan(H2fit));

a0 = [1 7];
[a] = lsqcurvefit(@myfun,a0,xdata,H2fit); %least square

%do the inverse of the function to get a predicted map of phi
Phipred=real(acos((datastack.H-a(2))/a(1))*90/pi());
datastack.Phipred=Phipred;

XPCcontourf(Phipred, 'title', "Simulated Declination Angle map",'cunits', "Declination angle /^{o}", ...
    'saveq',1)

XPCcontourf(datastack.Phirefl*180/pi()-datastack.Phipred, 'title',...
    "First correction Phi - Predicted Phi",'cunits', "Declination angle /^{o}", ...
    'saveq',1)

%% find the edges

edgeDetectorx = [1 2 1; 0 0 0; -1 -2 -1]; % Sobel edge detector.
edgeDetectory = [-1 0 1; -2 0 2; -1 0 1]; % Sobel edge detector.

Phipred2=padarray(datastack.Phipred,[4 4],'replicate'); % pad to remove edge effects
Phirefl2=padarray(datastack.Phirefl,[4 4],'replicate'); % pad to remove edge effects

Phiprededgesx = conv2(Phipred2*pi()/180, edgeDetectorx, 'same');
Phiedgesx = conv2(Phirefl2, edgeDetectorx, 'same');

Phiprededgesy = conv2(Phipred2*pi()/180, edgeDetectory, 'same');
Phiedgesy = conv2(Phirefl2, edgeDetectory, 'same');

Phiprededges=(180/pi)*(Phiprededgesx.^2+Phiprededgesy.^2).^0.5; %sum square
Phiprededges=Phiprededges(5:end-4,5:end-4); %remove pad

Phiedges=(180/pi)*(Phiedgesx.^2+Phiedgesy.^2).^0.5;
Phiedges=Phiedges(5:end-4,5:end-4); %remove pad


XPCcontourf(Phiprededges,'title',"Predicted Declination Angle Edges",...
    'cunits',"Gradient in declination angle, ^{o}/\mum",'saveq',1)

XPCcontourf(Phiedges,'title',"First correction Declination Angle Edges",...
    'cunits',"Gradient in declination angle, ^{o}/\mum",'saveq',1)

%% the second correction section 
transform = 'homography';%choose between homography and affine

if strcmp(transform,'affine')==1
    init=[eye(2) 0*ones(2,1)];%translation-only initialization
elseif strcmp(transform,'homography')==1
    init=eye(3);
    init(1:2,3) = 0;%translation-only initialization
end

NoI = 100; % number of iterations
NoL = 3;  % number of pyramid-levels
template_demo=Phiprededges; %the template is taken as havine GB in the right place
im_demo=padarray(Phiedges,[4 4],'replicate'); % the one to be "fixed"
%NOTE - padarray is necessary in order to make sure that any deviations
%don't cause a row/column of 0s on the edges due to lack of data. 
% This function does all the work
[~, final_warp, ~]=ecc(im_demo, template_demo, NoL, NoI, transform, init);
nx = 1:size(template_demo,2); %what should be the size of the final array:
ny = 1:size(template_demo,1);
% Using spatial_interp from the ECC functions is neat
Phiseccorr   = spatial_interp(double(padarray(datastack.Phi,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
phi1seccorr  = spatial_interp(double(padarray(datastack.phi1,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
phi2seccorr  = spatial_interp(double(padarray(datastack.phi2,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
phaseseccorr = spatial_interp(double(padarray(datastack.phase,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
BCebsdseccorr=  spatial_interp(double(padarray(datastack.BCebsd,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);

datastack.Phiseccorr   =Phiseccorr;%*180/pi();
datastack.phi1seccorr  =phi1seccorr;%*180/pi();
datastack.phi2seccorr  =phi2seccorr;%*180/pi();
datastack.phaseseccorr =phaseseccorr;
datastack.BCebsdseccorr=BCebsdseccorr;

datastack.Phireflseccorr=datastack.Phiseccorr;
datastack.Phireflseccorr(datastack.Phireflseccorr>(pi/2))=pi-datastack.Phireflseccorr(datastack.Phireflseccorr>(pi/2));
 
g_writeEBSDdata_seccorr(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'secondcorrected' currdate]),datastack)
 
%% plotting
 
XPCcontourf(datastack.Phiseccorr.*180/pi,'title',"Second Corrected Phi",...
    'cunits',"Declination angle /^{o}", 'saveq',1)
 
XPCcontourf(datastack.Phireflseccorr*180/pi()-datastack.Phipred,...
    'title',"Second Corrected Phi - Predicted Phi", 'cunits',"Declination angle /^{o}", 'saveq',1)
 
XPCcontourf(datastack.Phireflseccorr*180/pi()-datastack.Phirefl*180/pi(),...
    'title',"Second Corrected Phi - First Corrected Phi", 'cunits',"Declination angle /^{o}", 'saveq',1)
 
figure;
scatter(datastack.Phireflseccorr(:)*180/pi(),datastack.H(:),'x')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(datastack.H(:))-5*nanstd(datastack.H(:)) nanmean(datastack.H(:))+5*nanstd(datastack.H(:))])
xlim([0 90])
title('Declination angle against measured hardness')
figname=['SECCORR phi vs H_useedge' num2str(useedgeq) ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
 
 
%% Grain boundary analysis
disp('Running grain boundary analysis on second correction')
 
[datastack,ebsdCorrectedLGG,Grains]=g_dist2grainb_secondcorr(resultsdir, ebsdname,datastack,currdate,ebsd);
 


function F = myfun(a,data)
    x = data(:);
    F=a(1)*cos((pi()/90)*x)+a(2);
end

end
