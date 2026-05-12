function A_grad = grad_fun(A_mat, D)

%% This function computes the gradient in matrix form

m = size(A_mat, 1);
d = size(A_mat, 3);

A_grad = {};

for i = 1:m

    Asum = zeros(d, 1);
    normfacti = D(i);

    for j = 1:m

        atemp = squeeze(A_mat(i, j, :));

        Asum = Asum + atemp;

        A_grad{i, j} = (1 / normfacti) * atemp;

    end

    A_grad{i, i} = -(1 / normfacti) * Asum;

end

A_grad = cell2mat(A_grad);

end
