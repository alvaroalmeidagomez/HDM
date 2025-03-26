function [Tangvect  d ] = localPCA(tol,K,Mvector)

%%%%%%This function computes the local PCA with a certain tolerance tol
%%%%K is the number of points in the local neighbourhood
%%%%Indk

stemp=size(Mvector);
dimMan=[];
CeldSVD={};
TangOrg=[];
OrgIndx=[];


for i=1:stemp(1)
    Mtemp=Mvector(i,:,:);
    Mtemp=reshape(Mtemp,K,[]);

    [U,S,V] = svd(Mtemp);
    S=abs(diag(S));
    [S Indtemp]=sort(S,'descend');
    U=U(:,Indtemp);
    V=V(:,Indtemp);
 


    CeldSVD{end+1}={{U} {S} {V}}; %%A=USV


    dtemp=sum(S>tol);
    dimMan(end+1)=dtemp;
    
end

d=mode(dimMan);


Tangvect=zeros(stemp(1),d,stemp(3)); %%Array with tangvectors.


for i=1:stemp(1)

vtemp=CeldSVD{i};
vtemp=vtemp{3};
vtemp=cell2mat(vtemp);
vtemp=transpose(vtemp(:,1:d));
Tangvect(i,:,:)=vtemp;

end

  

end





