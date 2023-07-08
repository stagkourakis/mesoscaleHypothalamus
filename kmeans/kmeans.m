% MATLAB Script for Cluster Analysis and Visualization
%
% This script performs k-means clustering on a given dataset, evaluates
% the cluster assignment using silhouette analysis, recalculates the cluster 
% averages after removing clusters with a high correlation coefficient (CC), 
% and identifies correlations greater than a certain threshold (0.9 in this case). 
% Clusters identified as highly correlated are then merged. 
% 
% The script also visualizes the mean activity of all clusters and performs 
% hierarchical linkage clustering on the merged k-means clusters. The final
% clusters are visualized using dendrograms.
%
% Variables:
% - maxNumClusters: Maximum number of clusters to test in k-means clustering
% - explainedVarianceThreshold: Threshold for variance to be explained by clusters
% - repetitions: Number of times the k-means clustering will be performed
% - data: The input data to be clustered
% - clusterIndices: Indices of data for each cluster
% - clusterCentroids: Centroid of each cluster
% - clusterDistances: Sum of point-to-centroid distances in each cluster
% - numClusters: Number of clusters found by k-means
% - meanClusterActivity: Mean activity of each cluster
% - correlationMatrix: Matrix of correlation coefficients of mean cluster activities
% - pValueMatrix: Matrix of p-values corresponding to correlationMatrix
% - clusterIndicesMerged: clusterIndices after merging highly correlated clusters
% - numMergedClusters: Number of unique clusters after merging
% - clustersMerged: List of clusters that were merged due to high correlation
% - zScoreMergedClusters: z-score normalized merged clusters
% - Z: Linkage matrix resulting from hierarchical clustering
% - maxNumClusters: Maximum number of clusters to form in hierarchical clustering
% - T: Cluster assignments from hierarchical clustering
% - D: Pairwise distances of zScoreMergedClusters
%
% Please note that some of the called functions, like `kmeans_opt`, `plotMeanClusterActivity`, 
% `plotZscoreNeuralActivity`, and `plotDendrogram`, must be user-defined and available 
% in your MATLAB path for this script to run properly.

%%
%   Author: Stefanos Stagkourakis
%   California Institute of Technology
%   For questions please contact: stefanos.stagkourakis@gmail.com

%%

% Initial Parameters
maxNumClusters = 100;
explainedVarianceThreshold = 0.90;     % i.e., threshold for variance to be explained
repetitions = 10;

% Clear Variables
clear clusterIndices;
clear clusterCentroids;
clear clusterDistances;
clear numClusters;

% Perform k-means clustering on the data
[clusterIndices, clusterCentroids, clusterDistances, numClusters] = ...
    kmeans_opt(data, maxNumClusters, explainedVarianceThreshold, repetitions);

% Evaluate cluster assignment using silhouette method
[s, h] = silhouette(data, clusterIndices, 'Euclidean');

% Recalculate cluster averages after removing clusters with high correlation coefficient (CC)
meanClusterActivity = [];

for idx = 1:numClusters
    meanClusterActivity(idx,:) = mean(data(clusterIndices==idx,:));
end

% Calculate correlation of mean cluster activities
[correlationMatrix, pValueMatrix] = corr(meanClusterActivity');

% Find correlations higher than 0.9
[rhoDiagZHoriz1] = correlationMatrix - diag(diag(correlationMatrix));
[r,c] = find(triu(rhoDiagZHoriz1)>0.9);

% Merge the clusters identified with high correlation
clusterIndicesMerged = clusterIndices;
for idx = 1:numel(r)
    clusterIndicesMerged(clusterIndicesMerged==r(idx)) = c(idx);
end

% Calculate the number of unique clusters after merging
numMergedClusters = numel(unique(clusterIndicesMerged));

% Identify clusters that were merged due to high correlation
clustersMerged = [r c];

% Plot the mean activity of all clusters
plotMeanClusterActivity(clusterCentroids, 'S1 - RAW');

% Perform hierarchical linkage clustering on the k-means clusters
zScoreMergedClusters = zscore(numMergedClusters, [], 2);

% Plot z-scored neural activity
plotZscoreNeuralActivity(zScoreMergedClusters, 'S1, z-scored neural act');

% Perform linkage clustering
Z = linkage(zScoreMergedClusters(:,1:1624),'average','correlation');

% Find a maximum of N clusters in the data
maxNumClusters = 3;
T = cluster(Z,'maxclust', maxNumClusters);

% Create a dendrogram plot of Z. 
cutoff =0.91;

D = pdist(zScoreMergedClusters(:,1:1624));

plotDendrogram(Z, D, cutoff, maxNumClusters);
