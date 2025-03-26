function [Lap] = HodgeMatrix(KNeighpoints,MatrixA,tangv,d,korder,ep)


if(korder==d)
    if(d==1)
        Hb=computationMatrixH0(KNeighpoints,MatrixA,d);
    else
        Hb=computationMatrixH(KNeighpoints,MatrixA,tangv,d,korder-1);
    end
Lap=(korder)*Hb*transpose(Hb);
end

if(korder==0 && korder<d)
 H=computationMatrixH0(KNeighpoints,MatrixA,d);
 Lap=transpose(H)*H;   
end

if(korder==1 && korder<d)
 Hb=computationMatrixH0(KNeighpoints,MatrixA,d);
 Hf=computationMatrixH(KNeighpoints,MatrixA,tangv,d,korder);   
 Lap=((korder+1)*transpose(Hf)*Hf)+((korder)*Hb*transpose(Hb));
end


if(1<korder && korder<d)
Hb=computationMatrixH(KNeighpoints,MatrixA,tangv,d,korder-1);
Hf=computationMatrixH(KNeighpoints,MatrixA,tangv,d,korder);
Lap=((korder+1)*transpose(Hf)*Hf)+((korder)*Hb*transpose(Hb));
end

Lap=power(ep,-2)*Lap;

end

