function [sol] = embeddingfun(V,S,mtrun,d,korder,tm,Nsize)




Ncomb=nchoosek(d,korder);
Sn=diag(S);
Sn=power(Sn,tm/2);
Sn=Sn/max(Sn);
Sn=diag(Sn);


Vnew=V*Sn;
Vnew=Vnew(:,1:mtrun);

sol=zeros(Nsize,mtrun,mtrun);
for ind1=1:Nsize
   itemp1=(ind1-1)* Ncomb+1;
   itemp2=ind1*Ncomb;
  
   temM=Vnew(itemp1:itemp2,:);


   sol(ind1,:,:)=transpose(temM)*temM;
   
end



end

