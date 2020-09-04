%Data comes in as:
% 1     2       3       4       5       6
% x(Px)	y(Px)	phase   phi1    PHI     phi2

% Writes back clean data

function f_writeEBSDdata_seccorr(filename,datastack)


data=zeros(size(datastack.X,1)*size(datastack.Y,2),6);
data(:,1)=datastack.X(:);
data(:,2)=datastack.Y(:);
data(:,3)=datastack.phi1seccorr(:);
data(:,4)=datastack.Phiseccorr(:);
data(:,5)=datastack.phi2seccorr(:);
data(:,6)=datastack.phaseseccorr(:);

% Writes back clean data

% txt matrix
fileID = fopen(filename,'w');
formatSpec = '%f %f %f %f %f %f\n';


for idx = 1:size(data,1)
    
    % write txt matrix
    fprintf(fileID,formatSpec,data(idx,1),data(idx,2),data(idx,6),data(idx,3)*180/pi,data(idx,4)*180/pi,data(idx,5)*180/pi);
end

fclose(fileID);

end