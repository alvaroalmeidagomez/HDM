function [HM] = computationMatrixH(KNeighpoints,MatrixA,tangv,d,korder)

[KM,KPM] = orderPermutations(korder,d);

vtemp=size(MatrixA);
Npoints=vtemp(1);
K=vtemp(2);

Lk=nchoosek(d,korder);
Lkk=nchoosek(d,korder+1);


HM=zeros(Lkk*Npoints,Lk*Npoints);
for i=1:Npoints
    for j=1:K
        
      jback=j;
      j=KNeighpoints(j,i);
      Atemp=MatrixA{i,jback};
      
      Oi=tangv(i,:,:);
      Oi=reshape(Oi,d,[]);
      
      Oj=tangv(j,:,:);
      Oj=reshape(Oj,d,[]);
      Oj=transpose(Oj);

      Coefftemp=zeros(Lkk,Lk);  
  
   
 for L1=1:Lkk 
 for L2=1:Lk

     lkk=KPM(L1,:);
     lk=KM(L2,:);
     
    AL1=Atemp(lkk,1);
    Otemp=Oi(lkk,:)*Oj(:,lk);

    Coefftemp(L1,L2)=det([AL1 Otemp]);
 
 end
 end 
 ind1s=((i-1)*Lkk)+1;
 ind1e=i*Lkk;
 ind2s=((j-1)*Lk)+1;
 ind2e=j*Lk;
 
 HM(ind1s:ind1e,ind2s:ind2e)= Coefftemp;  
end
end


end

