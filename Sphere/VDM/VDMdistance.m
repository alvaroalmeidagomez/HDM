function [sol] = VDMdistance(EmbM,d,Nsize,ind1)


sol=zeros(Nsize,1);


    itemp1=(ind1-1)*d+1;
    itemp2=ind1*d;
   
   Mtemp1=EmbM(ind1,:,:);
   Mtemp1=reshape(Mtemp1,1,[]);
   
   for ind2=1:Nsize
   itemp3=(ind2-1)* d+1;
   itemp4=ind2*d;
   
   Mtemp2=EmbM(ind2,:,:);
   Mtemp2=reshape(Mtemp2,1,[]);
   
   var1=dot(Mtemp1,Mtemp1);
   var2=dot(Mtemp2,Mtemp2);
   var3=dot(Mtemp1,Mtemp2);
   
  
   
   sol(ind2)=sqrt(abs(real(var1+var2-2*var3)));
    end


end

