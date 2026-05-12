function Tang_M = tangvect(X)

    m = size(X,1);

    Tang_M = zeros(m,3,2);

    for i = 1:m

        vtemp = X(i,:);
        Tang_M(i,:,:) = null(vtemp);

    end

end