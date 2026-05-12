function [ED1] = ED1(Tang_M, A, D)

%%% This function computes the matrix of the exterior derivative

m = length(D);
l = size(A);
l = l(3);

ED11 = {};
ED21 = {};

for i = 1:m
    
    mtempi = squeeze(Tang_M(i,:,:));
    di = D(i);
    Asum = zeros(l,1);
    
    for j = 1:m
        
        mtempj = squeeze(Tang_M(j,:,:));
        
        A_temp = squeeze(A(i,j,:));
        Asum = Asum + A_temp;
        
        M_big_temp1 = [A_temp, mtempi' * mtempj(:,1)];
        M_big_temp2 = [A_temp, mtempi' * mtempj(:,2)];
        
        M_big = (1/di) * [det(M_big_temp1), det(M_big_temp2)];
        
        ED11{i,j} = M_big;
        ED21{i,j} = zeros(1,2);
        
    end
    
    M_big_temp1_n = [Asum, mtempi' * mtempi(:,1)];
    M_big_temp2_n = [Asum, mtempi' * mtempi(:,2)];
    
    M_big = (1/di) * [det(M_big_temp1_n), det(M_big_temp2_n)];
    
    ED21{i,i} = M_big;
    
end

ED11 = cell2mat(ED11);
ED21 = cell2mat(ED21);

ED1 = ED11 - ED21;

end

