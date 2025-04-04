clc
clear all
close all

u1=50; %%Number of sample points for the first coordinate of the Torus
u2=50; %% Number of sample points for the second coordinate of the Torus
K=30; %% k-near neighbourhood
tol=0.9; %% Tolerance to compute the rank in Local PCA
mtrun=3;
tm=1;
inpoint=1;

X=[];  %%X is the data set
Coor2d=[]; %%2d coordinates of the dataset
Vspace=linspace(0,1,u1); %%Definition of the vector with the first coordinate
Vspace2=linspace(0,1,u2); %%Definition of the vector with the second coordinate

for i1=1:u1
for i2=1:u2-1

   
utemp=2*pi*Vspace(i1); 
vtemp=pi*Vspace2(i2);
Coor2d(end+1,:)=[Vspace(i1), Vspace2(i2)];
X(end+1,:)=[cos(utemp)*sin(vtemp), sin(utemp)*sin(vtemp), cos(vtemp)];

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


%%%%%%%%%%%%% Load the data-set X %%%%%%%%%%%%%%%%

[KNeighpoints,Mvector,t] = CompMatrix(X,K); %%Calculating the matrix with K-nearest points

[tangv d ]=localPCA(tol,K,Mvector); %%LocaL PCA function

ConLap = vdm(tangv,X,t); %% Connection  Laplacian


[V,S]=svd(ConLap);

%%%%%%%%%%%Embedding Computation%%%%%%%%%%%%%

EmbM = embeddingfun(V,S,mtrun,d,tm,Nsize); %%Computing the truncated embedding


%%%%%%%%%%%Plotting the embedding%%%%%%%%%


for itempn1=1:mtrun
    for itempn2=itempn1:mtrun
        vtempn=EmbM(:,itempn1,itempn2);
        vtempn=reshape(vtempn,[],1);
        
        figure
        scatter(Coor2d(:,1),Coor2d(:,2),10,vtempn,'filled')
        title("The (" + itempn1 + "," + itempn2 + ") component of the vector diffusion map")
        xlabel('First coordinate')
        ylabel('Second coordinate')
        colorbar
        colormap jet
     
    end
end


%%%%%%%%%Hodge-Diffusion-Distance%%%%%%%%%

VDM=VDMdistance(EmbM,d,Nsize,inpoint);




figure
scatter3(X(:,1),X(:,2),X(:,3),10,VDM,'filled') %%Ploting the data-set X
hold on
scatter3(X(inpoint,1),X(inpoint,2),X(inpoint,3),100,'black','filled')
title('Vector diffusion maps distance')
colorbar
colormap jet


figure
scatter(Coor2d(:,1),Coor2d(:,2),10,VDM,'filled')
hold on
scatter(Coor2d(inpoint,1),Coor2d(inpoint,2),100,'black','filled')
title("Vector diffusion maps distance")
xlabel('First coordinate')
ylabel('Second coordinate')
colorbar
colormap jet


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
title("Vector diffusion maps")
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)


figure
scatter(DimRedM(:,1),DimRedM(:,2),10,colist,'filled')
colormap jet
xlabel('The (1,1)- coordinate')
ylabel('The (2,2)- coordinate')
title("Vector diffusion maps")
cb=colorbar;
ylabel(cb,'Sample point','FontSize',10)






