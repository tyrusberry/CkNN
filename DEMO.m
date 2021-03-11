clear;close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Generate data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);     %%% seed
numclus=3;  %%% number of clusteres
dim=2;      %%% dimension of the data set
gap=0.05;   %%% gap between the clusters

[X,clusterNums]=GenerateSpiralData(1000,dim,numclus,2,1,gap);



%%%%%%%%%%%%%%%%%%%% CkNN based clustering using CkNN %%%%%%%%%%%%%%%%%%%

[numclus,clusters,allMidpoints,allTransitions,ds,d] = CkNN(X);



%%%%%%%%%%%%%%%%%%%%% Plot the final clustering %%%%%%%%%%%%%%%%%%%%%%%%%

delta = ds(allTransitions(2))   %%% delta parameter that defines the CkNN graph construction
                                %%% extracted from the clustering output
figure(3);
plot(X(1,:),X(2,:),'.k','MarkerSize',10);
hold on;axis equal;axis off;set(gca,'FontSize',24);
T=size(X,2);
k=10;

for i=1:T
    for j=1:T
        %%% Apply the CkNN criterion to determine which edges to draw
        if (norm(X(:,i)-X(:,j)) <= delta*sqrt(d(k,i)*d(k,j)))
            if (clusters(i,2)==1)
                line([X(1,i) X(1,j)],[X(2,i) X(2,j)],'color','g');
            elseif (clusters(i,2)==2)
                line([X(1,i) X(1,j)],[X(2,i) X(2,j)],'color','c');
            else
                line([X(1,i) X(1,j)],[X(2,i) X(2,j)],'color','r');
            end
        end
    end
end
plot(X(1,:),X(2,:),'.k','MarkerSize',10);
title('Recovered Clusters','fontsize',20);