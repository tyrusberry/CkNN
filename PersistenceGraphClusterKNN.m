function [numClusters,clusterNumbers,allMidpoints,allTransitions,ds] = PersistenceGraphClusterKNN(ds,di,maxCluster,makeplots)

%%% Inputs
    %%% ds          - weighted distances to k nearest neighbors, N-by-k
    %%% di          - indices of the k nearest neighbors of each point, N-by-k
    %%% maxCluster  - maximum number of clusters to consider

%%% Outputs
    %%% numClus         -number of clusters corresponding to allMidpoints
    %%% clusterNumbers  - N-by-(maxCluster-1), for each number of clusters from 2 to maxCluster,
    %%%                     assigns a cluster number to each point
    %%% allMidpoints    -values of global scaling which give each # of clusters
    %%% allTransitions  -values of global scaling where # of clusters changes 
    %%% ds              -sorted list of weights of unique edges corresponding to edge 
    %%%                     numbers in allMidpoints and allTransitions
    
    if (nargin<4)
        makeplots=1;
    end
    N=size(ds,1);
    totEdges=N*(N-1)/2;
    k=size(ds,2);
    [ds,dii]= sort(ds(:));
    [is,~] = ind2sub([N k],dii);
    js = di(dii);
    alledges = [is js];
    alledges = sort(alledges,2);
    [~,ia] = unique(alledges,'rows');
    ia = sort(ia);
    alledges=alledges(ia,:);    %%% list of unique edges in order of weights
    ds = ds(ia);                %%% corresponding weights

    numEdges = size(alledges,1);
    numClus = 1;
    numClusSteps(1) = numClus;
    numEdgesSteps(1) = numEdges;
    
    count = 2;

    while (numClus < maxCluster)
        numEdges = floor(numEdges/2);
        g = graph(alledges(1:numEdges,1)',alledges(1:numEdges,2)');
        bins = conncomp(g);
        numClus = length(unique(bins));
        
        numClusSteps(count) = numClus;
        numEdgesSteps(count) = numEdges;
        
        count = count+1;
        if (makeplots)
            figure(makeplots);hold off;plot(100*numEdgesSteps/totEdges,numClusSteps,'ko');drawnow;
        end
    end

    minCluster = min(numClusSteps);
    
    %%% Find all the transition points via bisection method
    allTransitions = zeros(maxCluster-minCluster,1);
    for i=1:maxCluster-minCluster
        clusterNumber = i+minCluster;
        lb=max(numEdgesSteps(numClusSteps>=clusterNumber));
        ub=min(numEdgesSteps(numClusSteps<clusterNumber));
        allTransitions(i) = floor((lb+ub)/2);
        
        while (ub-lb > 1)

            numEdges = allTransitions(i);
            
            g = graph(alledges(1:numEdges,1)',alledges(1:numEdges,2)');
            bins = conncomp(g);
            numClus = length(unique(bins));

            numClusSteps(count) = numClus;
            numEdgesSteps(count) = numEdges;
            [numEdgesSteps,sinds]=sort(numEdgesSteps);
            numClusSteps = numClusSteps(sinds);
            count = count+1;
        
            if (numClus >= clusterNumber)
                lb = allTransitions(i);
            else
                ub = allTransitions(i);
            end
            
            %%% log scale bisection method
            allTransitions(i) = floor((lb+ub)/2);
            
            if (makeplots)
                figure(makeplots);hold off;plot(100*numEdgesSteps/totEdges,numClusSteps,'linewidth',2);
                hold on;plot(100*allTransitions(1:i)/totEdges,(1:i)+minCluster,'ro');
                xlim([min(100*numEdgesSteps/totEdges) max(100*allTransitions/totEdges)*1.5]);
                drawnow;
            end

        end
 
    end
    
    
    steps = maxCluster-minCluster-1;
    clusterNumbers = zeros(N,steps);
    numClusters = zeros(steps,1);
    allMidpoints = zeros(steps,1);
    
    if (minCluster > 1)
        allTransitions = [max(numEdgesSteps); allTransitions];
        steps = steps + 1;
    end
    
    %%% Cluster at the midpoints of the persistence regions
    for i=1:steps

        lb=allTransitions(i);
        ub=allTransitions(i+1);

        %%% log scale midpoint
        allMidpoints(i) = floor((lb+ub)/2);
        
        g = graph(alledges(1: allMidpoints(i),1)',alledges(1: allMidpoints(i),2)');
        bins = conncomp(g);
        clusterNumbers(:,i) = bins;
        numClusters(i) = length(unique(bins));

    end
    
    if (makeplots)
        percentTrans = 100*allTransitions/totEdges;
        set(gca,'Fontsize',18);
        xlabel('Percent of Edges','Fontsize',24);
        ylabel('Number of Clusters','Fontsize',24);
        ylim([0 maxCluster]);
        xlim([min(percentTrans) max(percentTrans)*1.2]);
        drawnow;   
    end