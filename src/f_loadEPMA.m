function [C, XC, YC] = f_loadEPMA(epmaname, filepath, epmabsename )
%Script for loading EPMA, as it was getting large. 

idx = strfind(epmaname,'.'); %get the extension and load ebsd

if strcmp(epmaname(idx+1:end),'png')
    [C,map] = imread(fullfile(filepath,epmaname),'png');
    if ~isempty(map)
        C = ind2rgb(C,map);
    end
    C = mean(C,3);
    C=im2double(C);
    C=rot90(rot90(C,1),2);
    xC = 1:size(C,1);
    yC = 1:size(C,2);
    XC=repmat(xC,[1,size(yC)]);
    XC=squeeze(XC);
    
    YC=repmat(yC,[1,size(xC)]);
    YC=squeeze(YC)';
    %smooth?
    C = smoothdata(C,2,'gaussian',2);
    C = smoothdata(C,1,'gaussian',2);    
elseif strcmp(epmaname(idx+1:end),'tiff') %if reading from a tiff, which has all the signal in it:
    C = imread(fullfile(filepath,epmaname),'tiff');
    try %if you've got a tiff maybe you've got the BSE image:
        BSEC = imread(fullfile(filepath,epmabsename),'tiff');
        C=cat(3,C, BSEC); %stack these two - the second layer is the BSE. this is a bit of a short sighted hack...
    catch
    end

    xC = 1:size(C,1);
    yC = 1:size(C,2);
    XC=repmat(xC,[1,size(yC)]);
    XC=squeeze(XC);
    
    YC=repmat(yC,[1,size(xC)]);
    YC=squeeze(YC)';
end

end
