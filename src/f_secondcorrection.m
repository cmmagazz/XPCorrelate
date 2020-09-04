%function [datastack]=f_secondcorrection(datastack,buffersize)
% This function performs a second correction to the Phi values. The way it
% works is by assuming that errors arise from grain boundaries being
% misplaced. The values near the grain boundaries are filtered, with
% distance set by BUFFERSIZE, and a fitting similar to lsqfitting is
% performed. The fit is then used to "predict" a phi for every point in the
% hardness map. This should not be nominally correct, but should have the
% grain boundaries in the right place. 

% This map is then used as a 
addpath(genpath('external'));

if exist('buffersize','var')~=1
    buffersize=input('What buffer size should be used?');
end

useedgeq=1;
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

figure;
contourf(datastack.X,datastack.Y,datastack.Phirefl*180/pi()-datastack.Phipred,45,'LineColor','None');
title('First correction Phi - Predicted Phi')
caxis([-10 10])
axis image
c=colorbar;
c.Label.String = 'Phi /deg';
figname=['Phi-Phipred' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')

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
if useedgeq==0
    template_demo=datastack.Phipred*pi()/180; %the template is taken as havine GB in the right place
    im_demo=padarray(datastack.Phirefl,[4 4],'replicate'); % the one to be "fixed"
    %NOTE - padarray is necessary in order to make sure that any deviations
    %don't cause a row/column of 0s on the edges due to lack of data. 
    % This function does all the work
    [~, final_warp, ~]=ecc(im_demo, template_demo, NoL, NoI, transform, init);
    nx = 1:size(template_demo,2); %what should be the size of the final array:
    ny = 1:size(template_demo,1);
    % Using spatial_interp from the ECC functions is neat
    Phiseccorr   = spatial_interp(double(im_demo), final_warp, 'nearest', transform, nx, ny);
    phi1seccorr  = spatial_interp(double(padarray(datastack.phi1,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
    phi2seccorr  = spatial_interp(double(padarray(datastack.phi2,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
    phaseseccorr = spatial_interp(double(padarray(datastack.phase,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
    BCebsdseccorr=  spatial_interp(double(padarray(datastack.BCebsd,[4 4],'replicate')), final_warp, 'nearest', transform, nx, ny);
elseif useedgeq==1
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
end




datastack.Phiseccorr   =Phiseccorr;%*180/pi();
datastack.phi1seccorr  =phi1seccorr;%*180/pi();
datastack.phi2seccorr  =phi2seccorr;%*180/pi();
datastack.phaseseccorr =phaseseccorr;
datastack.BCebsdseccorr=BCebsdseccorr;

datastack.Phireflseccorr=datastack.Phiseccorr;
datastack.Phireflseccorr(datastack.Phireflseccorr>(pi/2))=pi-datastack.Phireflseccorr(datastack.Phireflseccorr>(pi/2));

f_writeEBSDdata_seccorr(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'secondcorrected' currdate]),datastack)
[datastack,ebsdCorrectedLGG,Grains]=f_dist2grainb_secondcorr(resultsdir, ebsdname,datastack,currdate,ebsd);

%% plotting

figure;
contourf(datastack.X,datastack.Y,datastack.Phiseccorr,45,'LineColor','None');
title('Second Corrected Phi')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Phi /deg';
figname=['Second Corrected Phi_useedge' num2str(useedgeq) ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')


figure;
contourf(datastack.X,datastack.Y,datastack.Phireflseccorr*180/pi()-datastack.Phipred,45,'LineColor','None');
title('Second Corrected Phi - Predicted Ph')
xlabel('\mum')
ylabel('\mum')
axis image
caxis([-10 10])
c=colorbar;
c.Label.String = 'Phi /deg';
figname=['Phiseccorr-Phipred_useedge' num2str(useedgeq) ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')


figure;
scatter(datastack.Phireflseccorr(:)*180/pi(),datastack.H(:),'x')
xlabel('Declination angle /^{o}')
ylabel('Hardness /GPa')
ylim([nanmean(datastack.H(:))-5*nanstd(datastack.H(:)) nanmean(datastack.H(:))+5*nanstd(datastack.H(:))])
xlim([0 90])
title('Declination angle against measured hardness')
figname=['SECCORR phi vs H_useedge' num2str(useedgeq) ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')





function F = myfun(a,data)
    x = data(:);
    F=a(1)*cos((pi()/90)*x)+a(2);
end

%end
