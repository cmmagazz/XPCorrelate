function [x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,x_transREF,y_transREF]=f_ebsdwrongway(x_ebsd,y_ebsd,phi1,Phi,phi2,phase,BCebsd,EBSDREFq);
%If the ebsd was the wrong way round, pick points by hand to point it the
%right way up. 
%CMM 2020
ebsdfigww=figure;
if EBSDREFq==0 
    hplot=contourf(x_ebsd,y_ebsd,Phi,10,'LineColor','None');
elseif EBSDREFq==1
    hplot=contourf(x_ebsd,y_ebsd,phase,10,'LineColor','None');
elseif EBSDREFq==2
    hplot=contourf(x_ebsd,y_ebsd,BCebsd,10,'LineColor','None');
    caxis([nanmean(BCebsd(:))-nanstd(BCebsd(:)) nanmean(BCebsd(:))+nanstd(BCebsd(:))]) 
end
hold on
axis image
hold off
title({'Click on three points to determine direction:','first a point near (0,0), then (1,0), then (0,1)'})
[x_ww,y_ww] = getpts;
close(ebsdfigww)

%The way we'll fix this is by seeing how the points look in size. Calculate
%the delta between (0,0), (1,0), and (0,1).
dx_ww=[x_ww(2)-x_ww(1),x_ww(3)-x_ww(1)];
dy_ww=[y_ww(2)-y_ww(1),y_ww(3)-y_ww(1)];
%finding the highest value to see if the biggest one is negative or not.
[~,posx]=max(dx_ww);
[~,posy]=max(dy_ww);
%treat these as basis vectors: (this is a hack, if anyone knows a better way...
%there must be one!)

%The biggest value determines where the x/y direction points, and the 
if abs(dx_ww(1))> abs(dx_ww(2)) && abs(dy_ww(2))>abs(dy_ww(1)) && posx==1 && posy==2
    %if everythings fine
    %do nothing
elseif abs(dx_ww(1))> abs(dx_ww(2)) && abs(dy_ww(2))>abs(dy_ww(1)) && posx==2 && posy==2
    phi1=flipud(phi1);
    Phi=flipud(Phi);
    phi2=flipud(phi2);
    phase=flipud(phase);
    BCebsd=flipud(BCebsd);
elseif abs(dx_ww(1))> abs(dx_ww(2)) && abs(dy_ww(2))>abs(dy_ww(1)) && posx==1 && posy==1
    phi1=fliplr(phi1);
    Phi=fliplr(Phi);
    phi2=fliplr(phi2);
    phase=fliplr(phase);
    BCebsd=fliplr(BCebsd);
elseif abs(dx_ww(1))> abs(dx_ww(2)) && abs(dy_ww(2))>abs(dy_ww(1)) && posx==2 && posy==1
    phi1=flipud(fliplr(phi1));
    Phi=flipud(fliplr(Phi));
    phi2=flipud(fliplr(phi2));
    phase=flipud(fliplr(phase));
    BCebsd=flipud(fliplr(BCebsd));
    %these were the easy cases with only flips along x or y.
    %now we need flip along x=y
elseif abs(dx_ww(2))> abs(dx_ww(1)) && abs(dy_ww(1))>abs(dy_ww(2)) && posx==2 && posy==1
    x_ebsd=rot90(x_ebsd,1);
    y_ebsd=rot90(y_ebsd,1);
    Phi=rot90(Phi,1);
    phi1=rot90(phi1,1);
    phi2=rot90(phi2,1);
    phase=rot90(phase,1);
    BCebsd=rot90(BCebsd,1);
    x_ebsdt=y_ebsd; y_ebsdt=x_ebsd;
    x_ebsd=x_ebsdt; clearvars x_ebsdt; y_ebsd=y_ebsdt; clearvars y_ebsdt;
elseif abs(dx_ww(2))> abs(dx_ww(1)) && abs(dy_ww(1))>abs(dy_ww(2)) && posx==1 && posy==1
    phi1=flipud(phi1);
    Phi=flipud(Phi);
    phi2=flipud(phi2);
    phase=flipud(phase);
    BCebsd=flipud(BCebsd);
    x_ebsd=rot90(x_ebsd,1);
    y_ebsd=rot90(y_ebsd,1);
    Phi=rot90(Phi,1);
    phi1=rot90(phi1,1);
    phi2=rot90(phi2,1);
    phase=rot90(phase,1);
    BCebsd=rot90(BCebsd,1);
    x_ebsdt=y_ebsd; y_ebsdt=x_ebsd;
    x_ebsd=x_ebsdt; clearvars x_ebsdt; y_ebsd=y_ebsdt; clearvars y_ebsdt;
elseif abs(dx_ww(2))> abs(dx_ww(1)) && abs(dy_ww(1))>abs(dy_ww(2)) && posx==2 && posy==2
    phi1=fliplr(phi1);
    Phi=fliplr(Phi);
    phi2=fliplr(phi2);
    phase=fliplr(phase);
    BCebsd=fliplr(BCebsd);
    x_ebsd=rot90(x_ebsd,1);
    y_ebsd=rot90(y_ebsd,1);
    Phi=rot90(Phi,1);
    phi1=rot90(phi1,1);
    phi2=rot90(phi2,1);
    phase=rot90(phase,1);
    BCebsd=rot90(BCebsd,1);
    x_ebsdt=y_ebsd; y_ebsdt=x_ebsd;
    x_ebsd=x_ebsdt; clearvars x_ebsdt; y_ebsd=y_ebsdt; clearvars y_ebsdt;
elseif abs(dx_ww(2))> abs(dx_ww(1)) && abs(dy_ww(1))>abs(dy_ww(2)) && posx==1 && posy==2
    phi1=flipud(fliplr(phi1));
    Phi=flipud(fliplr(Phi));
    phi2=flipud(fliplr(phi2));
    phase=flipud(fliplr(phase));
    BCebsd=flipud(fliplr(BCebsd));
    x_ebsd=rot90(x_ebsd,1);
    y_ebsd=rot90(y_ebsd,1);
    Phi=rot90(Phi,1);
    phi1=rot90(phi1,1);
    phi2=rot90(phi2,1);
    phase=rot90(phase,1);
    BCebsd=rot90(BCebsd,1);
    x_ebsdt=y_ebsd; y_ebsdt=x_ebsd;
    x_ebsd=x_ebsdt; clearvars x_ebsdt; y_ebsd=y_ebsdt; clearvars y_ebsdt;
end

%correct for our standards of x and y
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
%re-centre on 0.
if x_ebsd(1,1)==0 || y_ebsd(1,1)==0 
    %do nothing
else
    x_ebsd= x_ebsd-x_ebsd(1,1);
    y_ebsd= y_ebsd-y_ebsd(1,1);
end
if x_ebsd(1,1)==0 && y_ebsd(1,1)==0 && x_ebsd(2,1)>x_ebsd(1,1) && y_ebsd(1,2)>y_ebsd(1,1)
    disp('ebsd zerod correctly')
else
    error('ebsd data not correctly structured')
end
%re-pick points for ebsd:
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
end