function [A_matrix, D] = compA(DS, Tang_M, t)

%%% This function computes the matrix A and the normalization vector D

stemp = size(DS);
stemp2 = size(Tang_M);

m = stemp(1);
d_dim_ext = stemp(2);
d_dim_int = stemp2(3);

A_matrix = zeros(m, m, d_dim_int);
D = zeros(m, 1);

for i = 1:m
    xi = DS(i,:);
    vtemp3 = squeeze(Tang_M(i,:,:));
    
    sum_acu = 0;
    
    for j = 1:m
        xj = DS(j,:);
        
        diff_vec = xj - xi;
        diff_vec = diff_vec';
        
        kernel_val = exp(-norm(diff_vec)^2 / (2 * t * t));
        
        tangent_proj = vtemp3' * diff_vec;
        weighted_val = kernel_val * tangent_proj;
        
        A_matrix(i, j, :) = weighted_val;
        sum_acu = sum_acu + kernel_val;
    end
    
    D(i) = sum_acu;
end

end

