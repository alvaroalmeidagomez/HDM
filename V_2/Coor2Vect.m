function Vecttang = Coor2Vect(Mcoord, Tang_M)

%% This function transforms 2D coordinates in each basis into a vector field (3D vectors)

% Reshape coordinate matrix into 2 rows
Mcoord2 = reshape(Mcoord, 2, []);

% Dimensions
m = size(Tang_M, 1);
n = size(Tang_M, 2);

% Output initialization
Vecttang = zeros(m, n);

for i = 1:m

    % Extract i-th tangent matrix (n x 2 or 2 x n depending on input)
    Mtemp = squeeze(Tang_M(i,:,:));

    % Apply transformation
    Vecttang(i,:) = Mtemp * Mcoord2(:, i);

end

end