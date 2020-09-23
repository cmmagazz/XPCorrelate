%Input deck for multi-dimensional spatially resolved analysis of EBSD,
%Express, and other maps. CMM 2019. 
%% 1
clear all
close all
clc
home
%adapted from MTEX - load mtex
try
  fid = fopen('VERSION','r');
  MTEXversion = fgetl(fid);
  fclose(fid);
  fprintf([' ' MTEXversion '  ']);
catch
    addpath Z:\CM\mtex-5.0.3
    startup_mtex
end

addpath(genpath('src'))
addpath(genpath('external'))

%% 2 User Inputs

filepath    = 'Z:\CM\18_OctEXPRESS\200727_CpTiGB\';  %Location of files
expressname = 'expresscptigb_5x5_+-82umspacing_40x40_2umspacing_3mNLC_Express_results.mat';
ebsdname    = 'preexpressmap.ctf'; %name of the ebsd with extension

epmaq       =  0; %epma analysis?
epmaname    = 'Map 1_ds1_O  Ka_O  Ka (Sp 1)_item2.tiff'; %name of epma file if needed 
epmabsename = 'Map 1_ds1Vs1 BSE Z_BSE Z_item1'; %name of epma BSE file if needed 

gdcalcq     = 1; %grain boundary analysis?
microscope  ='merlin2';%what microscope was the ebsd map taken on? 'evo', 'merlin', 'xbeam', or 'tescan'
primphase   ='Titanium';%what's the primary phase in the ebsd map? default uses 1st material phase
%ebsd analysis questions
EBSDREFq    = 0;%0 if just one phase: uses IPF for EBSD reference selection
                %1 if multiphase: uses phase map
                %2 if you want to use the pattern quality map
EPMArefq    = 1;%0 uses the chemical map for referencing
                %1 if present, uses the BSE map for referencing

hexmat      = 1; %if it's hexagonal, make reflect Phi values above 90 about 90deg
saveasebsdq = 1; %activate the writedata to a file that can be read by mtex
saveasfigq  = 0; %save images as .fig as well as .png?
%% 3 Loading and some manipulation
%network drive determination
NDD=filepath(1);
idx = strfind(ebsdname,'.'); %get the extension and load ebsd
if strcmp(ebsdname(idx+1:end),'h5')
    ebsd=loadEBSD_h5(fullfile(filepath,ebsdname));
else
    ebsd=loadEBSD(fullfile(filepath,ebsdname),ebsdname(idx+1:end));
end
filepathnew=filepath;
load(fullfile(filepath,expressname));%load the express data
filepath=filepathnew; %keep the user inputted filepath
%make sure the drive name is the same:
if strcmp(NDD,resultsdir(1)) && strcmp(NDD,filepath(1)) %take the resultsdir from the express input code
    %do nothing
else
    resultsdir=[NDD resultsdir(2:end)];
    filepath=[NDD filepath(2:end)];
end

currdate=datestr(datetime);
currdate=currdate(1:11);

resultsdirold=resultsdir;
idxres = strfind(resultsdirold,'\');
if idxres(end)==size(resultsdirold,2) %if there's a "\" at the end, dont take it
    idxres=idxres(end-2:end-1);
else
    idxres(1)=idxres(end); %if there isn't, take the last folder name
    idxres(2)=size(resultsdirold,2);
end
resultsdir=[filepath, resultsdirold(idxres(1)+1:idxres(2)),currdate];

%results directory creation 
if ~exist(resultsdir, 'dir')
   mkdir(resultsdir)
end

%Cropping if need be:
cropytop=0;
cropybot=0;
cropxright=0; 
cropxleft=0;
X=X(cropxleft+1:end-cropxright,cropytop+1:end-cropybot);
Y=Y(cropxleft+1:end-cropxright,cropytop+1:end-cropybot);
fullres=fullres(cropxleft+1:end-cropxright,cropytop+1:end-cropybot,:);


%% 4 Running the EBSD registration

[datastack]=f_fixEBSDdistortion(ebsd,X,Y,fullres,microscope,primphase, EBSDREFq,resultsdir);
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
% datastack.BCebsd          %Band Contrast (ebsd)

%CLEANING:
datastack.S(datastack.S>1e100)=0;     
datastack.D(datastack.D>1e100)=0;  
datastack.L(datastack.L>1e100)=0;    
datastack.M(datastack.M>1e100)=0;     
datastack.St(datastack.St>1e100)=0; 
datastack.H(datastack.H>1e100)=0;
%housekeeping:
X=datastack.X;
Y=datastack.Y;

if hexmat==1
    %first fix the phi to be between 0 and 90 FOR HCP?
    datastack.Phirefl=datastack.Phi;
    datastack.Phirefl(datastack.Phirefl>(pi/2))=pi-datastack.Phirefl(datastack.Phirefl>(pi/2));
end
%% 5 Saving figures
v_plotfig(datastack,resultsdir,ebsdname,saveasfigq)


if EBSDREFq==1
    v_multiphaseplotter(datastack,resultsdir,ebsdname)
end


%% 6 EPMA input
if epmaq==1
    cropq=0; %DO WE NEED TO CROP THE DATA? 1=yes 0=no
    %import the epma data (currently an image)
    [C, XC, YC] = f_loadEPMA(epmaname, filepath, epmabsename);
    
    %alignment
    [datastack]=f_fixEPMAdistortion(datastack,C,XC,YC,cropq, EPMArefq);
    
    %some smoothing
    datastack.EPMAOS = smoothdata(datastack.EPMAO,2,'gaussian',5.5);
    datastack.EPMAOS = smoothdata(datastack.EPMAOS,1,'gaussian',5.5);
    datastack.EPMAO = datastack.EPMAOS;
    
    v_plotfigepma(datastack,resultsdir, ebsdname)
end
%% 7 save as ebsd
if saveasebsdq==1
    f_writeEBSDdata(fullfile(resultsdir, [ebsdname(1:(max(size(ebsdname)-4))) 'corrected' currdate]),datastack)
end 

%% 8 grain distance analysis
if gdcalcq==1
    datastack=f_dist2grainb(resultsdir, ebsdname,datastack,currdate,ebsd);
end

%% 9 Save things
close all
save([fullfile(resultsdir,[filename(1:length(filename)-4) '_XPCorrelate_results' currdate]) '.mat']);

 
