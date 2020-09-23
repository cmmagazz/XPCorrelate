%% 2 sample t-test 
%a scrap script to get some stats if need be

figH=figure;
figH.Name='Modulus';
hplot=contourf(Xb,Yb,sampleb,45,'LineColor','None');
caxis([meanM-2*stdM meanM+2*stdM])
title('Final Modulus Map')
xlabel('\mum')
ylabel('\mum')
axis image
c=colorbar;

val1=1;
val2=50;
Xa=X(val1:val2,val1:val2);Ya=Y(val1:val2,val1:val2);
samplea=M(val1:val2,val1:val2);
Xb=X(end-val2:end,end-val2:end);Yb=Y(end-val2:end,end-val2:end);
sampleb=M(end-val2+1:end,end-val2+1:end);
[h,p,ci,stats] = ttest2(samplea(:),sampleb(:))


val3=62;
val4=61+50;
Xb=X(val3-50:val4-50,val3:val4);Yb=Y(val3-50:val4-50,val3:val4);
sampleb=M(val3-50:val4-50,val3:val4);
[h,p,ci,stats] = ttest2(samplea(:),sampleb(:))
hist=histogram(sampleb(:));
hold on
hist=histogram(samplea(:));

