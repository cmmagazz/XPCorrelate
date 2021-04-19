function v_multiphaseplotter(datastack,resultsdir,ebsdname)
%plots various things that are helpful if dealing with multiple phases and
%looking at phase-property relationships. CHECK THE VALUES FOR PHASE
%DISCRIMINATION 

    figphas=figure;
    figphas.Name='Phase';
    hplot=contourf(datastack.X,datastack.Y,datastack.phase,45,'LineColor','None');
    title('Final Phase Map')
    axis image
    c=colorbar;
    c.Label.String = 'Phase';
    figname=['Phase Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    
    figphasVH=figure;
    figphasVH.Name='phase vs H';
    scatter(datastack.phase(:),datastack.H(:))
    xlabel('Phase')
    ylabel('Hardness /GPa')
    %ylim([2 6])
    figname=['phase vs H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')

    % Find values of locations where the following conditions are satisfied
    where1 = datastack.phase >0.9 & datastack.phase < 1.1;
    where1p3 = datastack.phase > 1.2 & datastack.phase < 1.4;
    where1p7 = datastack.phase > 1.6 & datastack.phase < 1.8;
    where2 = datastack.phase > 1.9 & datastack.phase < 2.1;

    % Get means where H is in that phase
    meanH1 = nanmean(datastack.H(where1));
    meanH1p3 = nanmean(datastack.H(where1p3));
    meanH1p7 = nanmean(datastack.H(where1p7));
    meanH2 = nanmean(datastack.H(where2));
    
    stdevH1=nanstd(datastack.H(where1));
    stdevH1p3=nanstd(datastack.H(where1p3));
    stdevH1p7=nanstd(datastack.H(where1p7));
    stdevH2=nanstd(datastack.H(where2));
    
    H1=datastack.H(where1);
    H1p3=datastack.H(where1p3);
    H1p7=datastack.H(where1p7);
    H2=datastack.H(where2);
    
    Hexcl=datastack.H;
    Hexcl(where1)=NaN;
    H1=Hexcl; Hexcl=datastack.H;
    H1=Hexcl.*where1;
    H1(H1==0)=NaN;

    Hexcl(where1p3)=NaN;
    H1p3=Hexcl; Hexcl=datastack.H;

    Hexcl(where1p7)=NaN;
    H1p7=Hexcl; Hexcl=datastack.H;

    
    Hexcl(where2)=NaN;
    H2=Hexcl; Hexcl=datastack.H;
    H2=Hexcl.*where2;
    H2(H2==0)=NaN;

    M1=datastack.M(where1);
    M1p3=datastack.M(where1p3);
    M1p7=datastack.M(where1p7);
    M2=datastack.M(where2);

    phasevH=[1, meanH1,stdevH1; 2, meanH2,stdevH2];
    
    figphasVHM=figure;
    figphasVHM.Name='Phase vs Mean of H';
    scatter(phasevH(:,1),phasevH(:,2),'filled','o')
    hold on
    errorbar(phasevH(:,1),phasevH(:,2),phasevH(:,3),'LineStyle','none');
    xlabel('Phase')
    ylabel('Hardness /GPa')
    xlim([.9 2.1])
    xlabelnames = {'Phase 1';'Phase 2'};
    set(gca,'xtick',[1,2],'xticklabel',xlabelnames)
    hold off
    figname=['phase vs H Mean Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    saveas(gcf,fullfile(resultsdir, figname),'png')
    
%     %TO INCLUDE BITS IN BETWEEN PHASES WHEN IT'S INTERPOLATED
%     phasevH=[1, meanH1,stdevH1; 1.3, meanH1p3,stdevH1p3; 1.7, meanH1p7,stdevH1p7; 2, meanH2,stdevH2];
%     
%     figphasVHM=figure;
%     figphasVHM.Name='Phase vs Mean of H';
%     scatter(phasevH(:,1),phasevH(:,2),'filled','o')
%     hold on
%     errorbar(phasevH(:,1),phasevH(:,2),phasevH(:,3),'LineStyle','none');
%     xlabel('Phase')
%     ylabel('Hardness /GPa')
%     xlim([.9 2.1])
%     xlabelnames = {'Phase 1';'Phase 1 (near 2)';'Phase 2 (near 1)';'Phase 2'};
%     set(gca,'xtick',[1,1.3,1.7,2],'xticklabel',xlabelnames)
%     hold off
%     figname=['phase vs H Mean Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
%     saveas(gcf,fullfile(resultsdir, figname),'png')
%     
%     figphasMVH=figure;
%     figphasMVH.Name='Hardness vs Modulus';
%     scatter(M1(:),H1(:),15,'filld','r')
%     hold on
%     scatter(M1p3(:),H1p3(:),15,'filld','o','m')
%     scatter(M1p7(:),H1p7(:),15,'filld','o','g')
%     scatter(M2(:),H2(:),15,'filld','o','b')
%     hold off
%     xlabel('Modulus /GPa')
%     ylabel('Hardness /GPa')
%     legend('Phase 1','Phase 1 (near 2)','Phase 2 (near 1)','Phase 2','Location','Best');
%     %ylim([2 6])
%     figname=['Hardness vs Modulus Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
%     saveas(gcf,fullfile(resultsdir, figname),'png')
end