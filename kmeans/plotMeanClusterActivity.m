function plotMeanClusterActivity(meanClusterActivity, figureTitle)
    % Function to plot the mean activity of all clusters.

    figure('Name', figureTitle)
    imagesc(meanClusterActivity);
    colormap('jet');   % Change colormap to your preference.
    colorbar;
    xlabel('Behaviors (#)','FontSize',25,'FontWeight','bold')
    ylabel('Cluster ID (#)','FontSize',25,'FontWeight','bold')
    title('Active clusters vs behavior (not sorted)')
end
