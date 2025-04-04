function [sol] = embeddingfun(V,S,mtrun,d,tm,Nsize)




Sn=diag(S);
Sn=power(Sn,2*tm);
Sn=Sn/max(Sn);
Sn=diag(Sn);


Vnew=V*Sn;
Vnew=Vnew(:,1:mtrun);

sol=zeros(Nsize,mtrun,mtrun);
for ind1=1:Nsize
   itemp1=(ind1-1)*d+1;
   itemp2=ind1*d;
  
   temM=Vnew(itemp1:itemp2,:);


   sol(ind1,:,:)=transpose(temM)*temM;
   
end



end

