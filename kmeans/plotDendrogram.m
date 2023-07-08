function plotDendrogram(Z, D, cutoff, numClusters)
    % Function to plot a dendrogram.
    % Z: The linkage matrix.
    % D: The pairwise distances.
    % cutoff: The color threshold for the dendrogram.
    % numClusters: Maximum number of clusters to form.

    figure;
    subplot(3,1,1)
    dendrogram(Z,0,'ColorThreshold',cutoff);

    subplot(3,1,2)
    leafOrder = optimalleaforder(Z,D);
    dendrogram(Z,0,'reorder',leafOrder, 'ColorThreshold',cutoff);

    subplot(3,1,3)
    leafOrderSerial = 1:numClusters;
    dendrogram(Z,0,'reorder',leafOrderSerial, 'ColorThreshold',cutoff);
end
