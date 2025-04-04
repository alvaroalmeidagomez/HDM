function [ConLap] = vdm(tangv,X,ep)

stemp=size(tangv);
Nsample=stemp(1);
d=stemp(2);
dimext=stemp(3);

ConLap=zeros(Nsample*d,Nsample*d);

for ind1=1:Nsample
    A1temp=tangv(ind1,:,:);
    A1temp=reshape(A1temp,d,[]);
    itemp1=(ind1-1)*d+1;
    itemp2=ind1*d;
    for ind2=1:Nsample
    itemp3=(ind2-1)*d+1;
    itemp4=ind2*d;
        
    vtempNor=norm(X(ind1,:)-X(ind2,:));
    vtempNor=exp(-power(vtempNor,2)/(2*ep));
    
    A2temp=tangv(ind2,:,:);
    A2temp=reshape(A2temp,d,[]);
    
    Og=A1temp*transpose(A2temp);
    [Utemp,Stemp,Vtemp] = svd(Og);
    Og=Utemp*transpose(Vtemp);
    
    ConLap(itemp1:itemp2,itemp3:itemp4)=vtempNor*Og;
        
    end  
end

end

