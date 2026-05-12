function solfinal = VF(X, Tang_M)

%%% This function defines the vector field used in the heat equation

m = size(X, 1);

soltemp = zeros(size(X));

% Define base vector field in R^2
soltemp(:,1) = ones(m,1);
soltemp(:,2) =  ones(m,1);
soltemp(:,3) =  ones(m,1);

solfinal = {};

for i = 1:m

    % Vector field at point i
    vtemp = soltemp(i, :)';

    % Tangent matrix at point i
    Mtemp = squeeze(Tang_M(i,:,:));
    Mtemp = Mtemp';

    % Apply transformation
    solfinal{i,1} = Mtemp * vtemp;

end

solfinal = cell2mat(solfinal);

end