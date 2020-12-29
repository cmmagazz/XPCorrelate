function XPCcontourf(arraytoplot,varargin)
% A function to plot the way that I commonly use for XPcorrelate. 
% OBLIGATORY INPUTS: 
%       arraytoplot - an array the same size as X and Y to be plotted
%       
% OPTIONAL INPUTS
%       X - an X position array, where X(2,1) is greater than X(1,1)
%           Default is datastack.X
%       Y - a  Y position array, where Y(1,2) is greater than Y(1,1)
%           Default is datastack.Y
%       title - string to give figure a title IN DOUBLE QUOTES
%       limits - array to give limits in colour (e.g. [0 1])
%           Default is mean(arraytoplot) +/- 2 standard deviations
%       units - x/y units, IN DOUBLE QUOTES
%           Default is micron 
%       cunits - C units, placed near colourbar IN DOUBLE QUOTES
%       saveq - 1 if you want to save the figure
%           Default is no
%       resultsdir - default is global resultsdir, but can specify your own
%       IN DOUBLE QUOTES
%       saveasfigq - save as a fig? same as global by default
%
%   Example: XPCcontourf(datastack.H, 'title',"Hardness map", 'cunits',
%   "Hardness /GPa") 
%   This plots the H array in datastack, with the above title and colour
%   units, while assuming the rest as default.
% CMM 2020


try
    defaultx=evalin('base','datastack.X');
    defaulty=evalin('base','datastack.Y');
end
defaultfigtitle=string(inputname(1));
deflimits = [nanmean(arraytoplot(:))-2*nanstd(arraytoplot(:)) ...
    nanmean(arraytoplot(:))+2*nanstd(arraytoplot(:))];
defunits='\mum';
defcunits='';

defaultfigname=string(inputname(1));
defresultsdir = evalin('base','resultsdir');
defsaveq=0;
try
    defsaveasfigq = evalin('base','saveasfigq');
catch
    defsaveasfigq = 0;
end
p = inputParser;

validArrayPosNum = @(x) isnumeric(x);
addRequired(p,'arraytoplot',validArrayPosNum);

addOptional(p,'X',defaultx,validArrayPosNum);
addOptional(p,'Y',defaulty,validArrayPosNum);

addOptional(p,'limits',deflimits,validArrayPosNum);
addOptional(p,'title',defaultfigtitle,@isstring);
addOptional(p,'units',defunits,@isstring);
addOptional(p,'cunits',defcunits,@isstring);

addParameter(p,'saveq',defsaveq,@isscalar);
addOptional(p,'figname',defaultfigname,@isstring);
addOptional(p,'resultsdir',defresultsdir,@isstring);
addParameter(p,'saveasfigq',defsaveasfigq,@isscalar);

parse(p,arraytoplot,varargin{:});

figure;
hplot=contourf(p.Results.X,p.Results.Y,arraytoplot,45,'LineColor','None');
title(p.Results.title)
if size(unique(p.Results.limits),2)==2 %only set the Caxis if there's actually 2 values
    caxis(p.Results.limits)
end
xlabel(p.Results.units)
ylabel(p.Results.units)
axis image
c=colorbar;
if strlength(p.Results.cunits)>1 %only put c units if it's been inputted
    c.Label.String = p.Results.cunits;
end

if p.Results.saveq==1
    print(fullfile(p.Results.resultsdir, p.Results.figname),'-dpng','-r600')
    if p.Results.saveasfigq==1
        saveas(gcf,fullfile(p.Results.resultsdir, p.Results.figname),'fig')
    end
end

end