function [KNeighpoints,Mvector,ep] = CompMatrix(X,K)

%%% This function computes the K-nearest points

stemp=size(X);
DistM=zeros(stemp(1),stemp(1));
Mvector=zeros(stemp(1),K,stemp(2));
diagNorm=zeros(stemp(1),1);


for i1=1:stemp(1)
  for i2=1:stemp(1)
       vtemp=norm(X(i1,:)-X(i2,:));
       DistM(i1,i2)=vtemp;
 end
end

[B,IndKN] = sort(DistM);
KNeighpoints=IndKN(1:K,:);
ep=mean(B(2,:));

for i1=1:stemp(1)
    dtemp=0;
    sumtemp=zeros(1,stemp(2));
    for i2=1:K
        
        i2back=i2;
        i2=IndKN(i2,i1);
        vtemp=X(i2,:)-X(i1,:);
        vtemp1=norm(vtemp);
        exteVar=exp(-power(vtemp1,2)/(2*ep));
        Mvector(i1,i2back,:)=vtemp*exteVar;
        dtemp=dtemp+exteVar;
        sumtemp=sumtemp+vtemp*exteVar; 
    end
   
    Mvector(i1,1,:)=-sumtemp;
    diagNorm(i1)=dtemp;
   % Mvector(i1,:,:)=Mvector(i1,:,:)/dtemp;
 
    

end

for i1=1:stemp(1)
    for i2=1:K
        i2back=i2;
        i2=IndKN(i2,i1);
       Mvector(i1,i2back,:)=Mvector(i1,i2back,:)/sqrt( diagNorm(i1)*diagNorm(i2));
    end
end
 
end