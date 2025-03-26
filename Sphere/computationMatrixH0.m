function [HM] = computationMatrixH0(KNeighpoints,MatrixA,d)

vtemp=size(MatrixA);
Npoints=vtemp(1);
K=vtemp(2);

HM=zeros(d*Npoints,Npoints);
for i=1:Npoints
    for j=1:K
        
      jback=j;
      j=KNeighpoints(j,i);
      Atemp=MatrixA{i,jback};

 ind1s=(i-1)*d+1;
 ind1e=i*d;
 ind2s=(j-1)+1;
 ind2e=j;
 
 HM(ind1s:ind1e,ind2s:ind2e)= Atemp;  
end
end


end