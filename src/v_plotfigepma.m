function v_plotfigepma(datastack,resultsdir, ebsdname)
    resolution=['-r' num2str(600)];

    figFC=figure;
    figFC.Name='FC';
    hplot=contourf(datastack.X,datastack.Y,datastack.EPMAO,45,'LineColor','None');
    title('Registered EPMA map')
    %caxis([0 20])
    xlabel('\mum')
    ylabel('\mum')
    axis image
    c=colorbar;
    c.Label.String = 'Oxygen Concentration /arb units';
    figname=['FC Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    print(fullfile(resultsdir, figname),'-dpng',resolution)
    
    try
        figepbse=figure;
        figepbse.Name='EPMABSE';
        hplot=contourf(datastack.X,datastack.Y,datastack.EPMABSE,45,'LineColor','None');
        title('Final EPMA BSE Map')
        caxis([95 105])
        axis image
        c=colorbar;
        c.Label.String = 'arb units';
        figname=['EPMABSE Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
        print(fullfile(resultsdir, figname),'-dpng',resolution)
    end
    
    figFCVH=figure;
    figFCVH.Name='FC vs H';
    try
        scatter(datastack.EPMAO(:),datastack.Hmat(:))
    catch
        scatter(datastack.EPMAO(:),datastack.H(:))
    end
    xlabel('Oxygen Concentration /arb units')
    ylabel('Hardness /GPa')
    title({'EPMA obtained oxygen presence against';'measured nanoindentation hardness'})
    xlim([0 30])
    figname=['FC vs H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    print(fullfile(resultsdir, figname),'-dpng',resolution)

    figFCVVPH=figure;
    figFCVVPH.Name='FC vs phi vs H';
    try
        scatter3(datastack.EPMAO(:),datastack.Phi(:)*180/pi,datastack.Hmat(:))
    catch
        scatter3(datastack.EPMAO(:),datastack.Phi(:)*180/pi,datastack.H(:))
    end
    xlabel('Oxygen Concentration /arb units')
    ylabel('Declination angle /^{o}')
    title({'EBSD declination angle against EPMA oxygen concentration';' against measured nanoindentation hardness'})
    ylim([0 90])
    xlim([0 30])
    zlabel('Hardness /GPa')
    figname=['FC vs phi vs H Figure ' ebsdname(1:(max(size(ebsdname)-4)))];
    print(fullfile(resultsdir, figname),'-dpng',resolution)
    close all
end

