function v_multiphaseplotter(datastack,resultsdir,ebsdname)
%plots various things that are helpful if dealing with multiple phases and
%looking at phase-property relationships. CHECK THE VALUES FOR PHASE
%DISCRIMINATION 

    figphas=figure;
    figphas.Name='Phase';
    hplot=contourf(X,Y,datastack(:,:,12),45,'LineColor','None');
    title('Final Phase Map')
    axis image
    c=colorbar;
    c.Label.String = 'Phase';
    figname=['Phase Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    saveas(gcf,fullfile(resultsdir, figname),'fig')
    
    FphaseGr=datastack(:,:,12);
    figphasVH=figure;
    figphasVH.Name='phase vs H';
    scatter(FphaseGr(:),H(:))
    xlabel('Phase')
    ylabel('Hardness /GPa')
    %ylim([2 6])
    figname=['phase vs H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    saveas(gcf,fullfile(resultsdir, figname),'fig')    

    % Find values of locations where the following conditions are satisfied
    where1 = FphaseGr >0.9 & FphaseGr < 1.1;
    where1p3 = FphaseGr > 1.2 & FphaseGr < 1.4;
    where1p7 = FphaseGr > 1.6 & FphaseGr < 1.8;
    where2 = FphaseGr > 1.9 & FphaseGr < 2.1;

    % Get means where H is in that phase
    meanH1 = nanmean(H(where1));
    meanH1p3 = mean(H(where1p3));
    meanH1p7 = mean(H(where1p7));
    meanH2 = mean(H(where2));
    
    stdevH1=nanstd(H(where1));
    stdevH1p3=nanstd(H(where1p3));
    stdevH1p7=nanstd(H(where1p7));
    stdevH2=nanstd(H(where2));
    
    H1=H(where1);
    H1p3=H(where1p3);
    H1p7=H(where1p7);
    H2=H(where2);
    
    Hexcl=H;
    Hexcl(where1)=NaN;
    H1=Hexcl; Hexcl=H;
    H1=Hexcl.*where1;
    H1(H1==0)=NaN;

    Hexcl(where1p3)=NaN;
    H1p3=Hexcl; Hexcl=H;

    Hexcl(where1p7)=NaN;
    H1p7=Hexcl; Hexcl=H;

    
    Hexcl(where2)=NaN;
    H2=Hexcl; Hexcl=H;
    H2=Hexcl.*where2;
    H2(H2==0)=NaN;

    M1=M(where1);
    M1p3=M(where1p3);
    M1p7=M(where1p7);
    M2=M(where2);

    phasevH=[1, meanH1,stdevH1; 1.3, meanH1p3,stdevH1p3; 1.7, meanH1p7,stdevH1p7; 2, meanH2,stdevH2];
    
    figphasVHM=figure;
    figphasVHM.Name='Phase vs Mean of H';
    scatter(phasevH(:,1),phasevH(:,2),'filled','o')
    hold on
    errorbar(phasevH(:,1),phasevH(:,2),phasevH(:,3),'LineStyle','none');
    xlabel('Phase')
    ylabel('Hardness /GPa')
    xlim([.9 2.1])
    xlabelnames = {'Phase 1';'Phase 1 (near 2)';'Phase 2 (near 1)';'Phase 2'};
    set(gca,'xtick',[1,1.3,1.7,2],'xticklabel',xlabelnames)
    hold off
    figname=['phase vs H Mean Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    saveas(gcf,fullfile(resultsdir, figname),'fig')
    
    figphasMVH=figure;
    figphasMVH.Name='Hardness vs Modulus';
    scatter(M1(:),H1(:),15,'filld','r')
    hold on
    scatter(M1p3(:),H1p3(:),15,'filld','o','m')
    scatter(M1p7(:),H1p7(:),15,'filld','o','g')
    scatter(M2(:),H2(:),15,'filld','o','b')
    hold off
    xlabel('Modulus /GPa')
    ylabel('Hardness /GPa')
    legend('Phase 1','Phase 1 (near 2)','Phase 2 (near 1)','Phase 2','Location','Best');
    %ylim([2 6])
    figname=['Hardness vs Modulus Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    saveas(gcf,fullfile(resultsdir, figname),'fig')
end