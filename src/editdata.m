% File just for editting data that doesn't fit anywhere else
% CMM 2019


%% deleting some overlaying rows

for i =1:7
   H(:,i*41+1)=NaN;
end
H2=H;
H(33:end,:)=NaN;
H(:,1:25)=NaN;
H(:,245)=NaN;


figure;
hplot=contourf(X,Y,H,45,'LineColor','None');
%caxis([3 9.5]);
title('Hardness')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;
c.Label.String = 'Hardness (GPa)';



%% some ebsd playing
%{

ebsd1=loadEBSD(fullfile(filepath,'Ti834Map11.ctf'),'ctf');
ebsd2=loadEBSD(fullfile(filepath,'Ti834Map12.ctf'),'ctf');
ebsd3=loadEBSD(fullfile(filepath,'Ti834Map13.ctf'),'ctf');
ebsd4=loadEBSD(fullfile(filepath,'Ti834Map14.ctf'),'ctf');


dx12=-9.5;
dy12=-.5;
max_ebsd1 = max(ebsd1.x);
ebsd2s = shift(ebsd2,[max_ebsd1+dx12,+dy12]);
ebsd12=[ebsd1 ebsd2s];
plot(ebsd12,ebsd12.prop.bc); colormap('gray')
%
dx123=-3.5;
dy123=0;
max_ebsd12 = max(ebsd12.x);
ebsd3s = shift(ebsd3,[max_ebsd12+dx123,+dy123]);
ebsd123=[ebsd12 ebsd3s];
plot(ebsd123,ebsd123.prop.bc); colormap('gray')

%
dx1234=-8;
dy1234=-0.75;
max_ebsd123 = max(ebsd123.x);
ebsd4s = shift(ebsd4,[max_ebsd123+dx1234,+dy1234]);
ebsd1234=[ebsd123 ebsd4s];
plot(ebsd1234,ebsd1234.prop.bc); colormap('gray')

setMTEXpref('xAxisDirection','west');
setMTEXpref('zAxisDirection','outOfPlane');

region = [0 0 33.2 5.48]*10^1;
%plot(ebsd1234)
%rectangle('position',region,'edgecolor','r','linewidth',2)

condition = inpolygon(ebsd1234,region);
ebsdfinal = ebsd1234(condition);
plot(ebsdfinal,ebsdfinal.prop.bc); colormap('gray')
%}