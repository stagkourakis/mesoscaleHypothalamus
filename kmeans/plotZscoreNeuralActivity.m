function plotZscoreNeuralActivity(zScoreData, figureTitle)
    % Function to plot the z-scored neural activity.

    figure('Name', figureTitle)
    imagesc(zScoreData);
    colormap('jet');   % Change colormap to your preference.
    caxis([-2 2])
    colorbar;
    xlabel('Behaviors (#)','FontSize',25,'FontWeight','bold')
    ylabel('Cluster ID (#)','FontSize',25,'FontWeight','bold')
    title('Active clusters vs behavior (not sorted)')
end
