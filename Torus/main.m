clc
clear all
close all

%%%%%%%%%%Definition of the data-set (Torus)%%%%%%%%%%

u1=50; %%Number of sample points for the first coordinate of the Torus
u2=50; %% Number of sample points for the second coordinate of the Torus
K=30; %% k-near neighbourhood
tol=0.9; %% Tolerance to compute the rank in Local PCA
korder=2; %% Order (integer) to compute the k-th Hodge Laplacian
mtrun=3; %%Truncation order
tm=1;  %% Time of the Markov Chain


X=[];  %%X is the data set
Coor2d=[]; %%2d coordinates of the dataset
Vspace=linspace(-1/2,1/2,u1); %%Definition of the vector with the first coordinate
Vspace2=linspace(-1/2,1/2,u2); %%Definition of the vector with the second coordinate

for i1=1:u1-1
for i2=1:u2
 
utemp=2*pi*Vspace(i1);
vtemp=2*pi*Vspace2(i2);
Coor2d(end+1,:)=[Vspace(i1), Vspace2(i2)];
X(end+1,:)=[(2+cos(vtemp))*cos(utemp), (2+cos(vtemp))*sin(utemp), sin(vtemp)];

    end
end

[X,i1x,i2x]=unique(X,'stable','rows'); %% Removing duplicate sample points
Coor2d=Coor2d(i1x,:);   %% Computing the coordinates of the data-set
[C1,C2] = meshgrid(Coor2d(:,1),Coor2d(:,2));
stemp=size(X);
colist=1:stemp(1);
Nsize=stemp(1);

%%%%%%%%%%%%%%%%%Plot Data-set %%%%%%%%%%%%%%%%%%%

figure
scatter3(X(:,1),X(:,2),X(:,3),10,colist,'filled') %%Ploting the data-set X
title('Plot of the data-set X')
colormap jet
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)

figure
scatter(Coor2d(:,1),Coor2d(:,2),10,colist,'filled')
title('Plot of the parametrization system of the data-set X')
xlabel('First coordinate')
ylabel('Second coordinate')
colormap jet
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)

%%%%%%%%%%%%% Load the data-set X %%%%%%%%%%%%%%%%


tttttemp=tic;

[KNeighpoints,Mvector,t] = CompMatrix(X,K); %%Calculating the matrix with K-nearest points

[tangv d ]=localPCA(tol,K,Mvector); %%LocaL PCA function

MatrixA= compMatrixA(Mvector,tangv); %Auxiliar function to compute the Hodge Laplacian

LaplacianM=HodgeMatrix(KNeighpoints,MatrixA,tangv,d,korder,t); %% Function to compute the Hodge Laplacian

[V,S]=svd(LaplacianM); %% Computing the SVD decomposition of the Hodgee-Laplacian.


%%%%%%%%% End Hodge-Laplacian Computation %%%%%%%%%

%%%%%%%%%%%Embedding Computation%%%%%%%%%%%%%

EmbM = embeddingfun(V,S,mtrun,d,korder,tm,Nsize); %%Computing the truncated embedding


%%%%%%%%%%%Plotting the embedding%%%%%%%%%


for itempn1=1:mtrun
    for itempn2=itempn1:mtrun
        vtempn=EmbM(:,itempn1,itempn2);
        vtempn=reshape(vtempn,[],1);
        
        figure
        scatter(Coor2d(:,1),Coor2d(:,2),10,vtempn,'filled')
        title("The (" + itempn1 + "," + itempn2 + ") component of the Hodge diffusion map")
        xlabel('First coordinate')
        ylabel('Second coordinate')
        colorbar
        colormap jet

     
    end
end


%%%%%%%%Diagonal embedding
DimRedM=zeros(Nsize,mtrun);
for itempn1=1:mtrun
        vtempn=EmbM(:,itempn1,itempn1);
        DimRedM(:,itempn1)=reshape(vtempn,[],1);
end
        
figure
scatter3(DimRedM(:,1),DimRedM(:,2),DimRedM(:,3),10,colist,'filled')
colormap jet
xlabel('The (1,1)- coordinate')
ylabel('The (2,2)- coordinate')
zlabel('The (3,3)- coordinate')
title("Hodge diffusion maps")
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)


figure
scatter(DimRedM(:,1),DimRedM(:,2),10,colist,'filled')
colormap jet
xlabel('The (1,1)- coordinate')
ylabel('The (2,2)- coordinate')
title("Hodge diffusion maps")
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)



ttn=toc(tttttemp);

 ppri=['The dimension of the manifold is ',num2str(d), '.']; %%Print Algorithm progress
 disp(ppri)
 ppri=['Hodge diffusion maps completed in ',num2str(ttn), ' Seconds.']; %%Print Algorithm progress
 disp(ppri)
