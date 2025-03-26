function [A] = compMatrixA(Mvector,tangv)

stemp=size(tangv);
N=stemp(1);
d=stemp(2);
stemp=size(Mvector);
K=stemp(2);

A={};

for i=1:N
for j=1:K
  Oi=tangv(i,:,:);
  Oi=reshape(Oi,d,[]);
 
  Diftemp=Mvector(i,j,:);
  Diftemp=reshape(Diftemp,[],1);

  A{i,j}=Oi*Diftemp;
  
end
end

end

