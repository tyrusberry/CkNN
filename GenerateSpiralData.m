function [data,clusterNums]=GenerateSpiralData(N,dim,clusters,width,length,minrad)
%%% Generate data on spiral shaped clusters in dim dimensions

if (nargin<6) minrad=.4;    end
if (nargin<5) length=1;     end
if (nargin<4) width=1;      end
if (nargin<3) clusters=3;   end
if (nargin<2) dim=2;        end
if (nargin<1) N=5000;       end

ang=(0:(clusters-1))*pi/clusters;
if (dim==2) ang=2*ang; end

clusterNums = ceil(clusters*rand(N,1));

%%% allphase determines which spiral each point will be on
allphase = ang(clusterNums);

%%% allradii determines how far down the spiral the point will be
allradii = minrad+length*abs(randn(N,1))';

%%% x defines the centerlines of the spirals, dx defines the tangent vector
for i=1:dim-1
    x(i,:) = allradii.*cos(allradii+allphase).*(sin(allradii+allphase).^(i-1));
    dx(i,:) = cos(allradii+allphase).*(sin(allradii+allphase).^(i-1));
    dx(i,:) = dx(i,:) - allradii.*sin(allradii+allphase).*(sin(allradii+allphase).^(i-1));
    dx(i,:) = dx(i,:) + (i-1)*allradii.*(cos(allradii+allphase).^2).*(sin(allradii+allphase).^(i-2));
end
x(dim,:) =  allradii.*(sin(allradii+allphase).^(dim-1));
dx(dim,:) =  (sin(allradii+allphase).^(dim-1)) + (dim-1)*allradii.*(sin(allradii+allphase).^(dim-2)).*cos(allradii+allphase);
dx = dx./(repmat(sqrt(sum(dx.^2,1)),dim,1));  %%% unit tangent vector

%%% z are random points on the n-1 dimensional sphere
z = randn(dim,N);
z = z./(repmat(sqrt(sum(z.^2,1)),dim,1));
z = z - repmat(sum(z.*dx,1),dim,1).*dx;

%%% allwidth determines how far from the centerline each point is
allwidth = width*(pi/clusters/dim/2)*min(allradii,1).*(rand(1,N));
%allwidth = width*(pi/spokes/n/2)*log(1+allradii).*(rand(1,N));

%%% z are random points on n-1 ball orthogonal to dx, with radius
%%% proportional to width and distance from x to the origin (conical)
z = repmat(allwidth,dim,1).*z./(repmat(sqrt(sum(z.^2,1)),dim,1));

%%% the final data points is formed by placing the point z, centered at x,
%%% in the hyperplane orthogonal to dx
data = x + z;

if (dim==2)
    figure(2);scatter(data(1,:),data(2,:),15,clusterNums,'filled')
    mycolors = [1 0 0;0 1 0;0 1 1];
    title('True Clusters','fontsize',20);
    colormap(mycolors);
else  
    figure(1);scatter3(x(1,:),x(2,:),x(3,:),20,clusterNums,'filled')
    figure(2);scatter3(data(1,:),data(2,:),data(3,:),20,clusterNums,'filled')
    mycolors = [1 0 0;0 1 0;0 1 1];
    title('True Clusters','fontsize',20);
    colormap(mycolors);
end

end

 