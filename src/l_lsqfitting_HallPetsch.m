%% Function to run some least squares fitting to data on hall petch
%CMM 2019
%REQUIRES the binning in p_PCDM 

a0 = [max(Have(:)) -1 Have(end)];

xdata=sqrt(gbszave(:)/pi*4); %model the size as an effective radius
%turn into column of data and remove nan

[a, resnorm,residual,~,~,~,J] = lsqcurvefit(@myfun,a0,xdata,Have); %least square
hrsq=1-sum(residual.^2)/sum((Have-mean(Have)).^2); %1-sum squared regression/sum of squares

[ypred,delta] = nlpredci(@myfun,xdata,a,residual,'Jacobian',J,'PredOpt','observation'); %do a prediction
lower = ypred-delta;
upper = ypred+delta;
dataupperlower=[xdata,upper,lower];

figure;
data2=linspace(min(xdata),max(xdata));
F=a(1)*data2.^(a(2))+a(3);
plot(data2,F,'r','LineWidth',2)
hold on
dataupperlower=sortrows(dataupperlower);
plot(dataupperlower(:,1),[dataupperlower(:,2),dataupperlower(:,3)],'r:', 'linewidth',2)
scatter(xdata(:),Have(:),'bx')
txt = (['H = ',num2str(a(1),'%4.2f'), '*Grain Size ^{' ,num2str(a(2),'%4.2f'),'} + ', num2str(a(3),'%4.2f')]); %THIS SHIT IS WRONG
text(1, 6.7, txt);%, 'FontSize', 10, 'Color', 'k');
legend({'Exponential Fit','95% confidence interval'},'Location','northeast')
hold off 
xlabel('Grain Effective Radius /\mum')
ylabel('Hardness /GPa')
ylim([6.6 7.7])
%xlim([0 90])
title('Grain Size against Measured Nanoindenter Hardness')
figname=['GSize vs H FIT' ebsdname(1:(max(size(ebsdname)-4)))];
saveas(gcf,fullfile(resultsdir, figname),'png')
saveas(gcf,fullfile(resultsdir, figname),'fig')



%% functions

function F = myfun(a,data)
    x = data(:);
    F=a(1)*x.^(a(2))+a(3);
end