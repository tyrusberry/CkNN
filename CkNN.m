function [numClusters,clusterNumbers,allMidpoints,allTransitions,ds,d] = CkNN(X,maxClus,m,k,makePlot)

%%% Inputs
    %%% X           - n-by-N matrix representing N data points in R^n
    %%% maxCluster  - maximum number of clusters to consider
    %%% m           - number of nearest neighbors to compute
    %%% k           - number of nearest neighbors to use in CkNN normalization


%%% Outputs
    %%% numClus         -number of clusters corresponding to allMidpoints
    %%% clusterNumbers  - N-by-(maxCluster-1), for each number of clusters from 2 to maxCluster,
    %%%                     assigns a cluster number to each point
    %%% allMidpoints    -values of global scaling which give each # of clusters
    %%% allTransitions  -values of global scaling where # of clusters changes 
    %%% ds              -sorted list of weights of unique edges corresponding to edge 
    %%%                     numbers in allMidpoints and allTransitions
    %%% d               -distances to the m nearest neighbors
    
    N=size(X,2);
    
    if (nargin<5) makePlot=1;  end
    if (nargin<4) maxClus = 10; end
    if (nargin<3) m = min(500,N); end
    if (nargin<2) k = 10; end

    [d,inds] = pdist2(X',X','euclidean','smallest',m); %%% Compute the m-NN of each point

    qt = d(k,:);                                %%% Compute the distance to the k-th NN
    CkNN = d./sqrt(repmat(qt,m,1).*qt(inds));   %%% This is the CkNN normalization

    %%%%%%%%%%%%%%%%% Persistence based clustering using CkNN %%%%%%%%%%%%%%%%

    [numClusters,clusterNumbers,allMidpoints,allTransitions,ds] = PersistenceGraphClusterKNN(CkNN',inds',maxClus,makePlot);


end

